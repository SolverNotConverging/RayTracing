#include "PathRefinement.hpp"

#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace {
bool positive_finite(double value) { return std::isfinite(value) && value > 0.0; }

bool usable(const SequenceEvaluation &geometry) {
    return geometry.status_ == SequenceStatus::Valid && geometry.receiver_ &&
        geometry.receiver_->residual_.allFinite() &&
        std::isfinite(geometry.receiver_->missDistance_) &&
        std::isfinite(geometry.receiver_->pathDistance_);
}

bool collapsing_corner(const SequenceEvaluation &geometry, double tolerance) {
    for (std::size_t i = 1; i < geometry.reflections_.size(); ++i) {
        const auto &a = geometry.reflections_[i - 1];
        const auto &b = geometry.reflections_[i];
        if ((a.position_ - b.position_).norm() <= tolerance &&
            std::min((a.normal_ - b.normal_).norm(), (a.normal_ + b.normal_).norm()) > equivalentNormalTolerance)
            return true;
    }
    return false;
}

bool same_path(const RefinementResult &a, const RefinementResult &b, const DeduplicationOptions &options) {
    const auto &ga = a.geometry_;
    const auto &gb = b.geometry_;
    if (ga.reflections_.size() != gb.reflections_.size()) return false;
    if ((a.transmitterPosition_ - b.transmitterPosition_).norm() > options.positionTolerance_ ||
        (a.receiverPosition_ - b.receiverPosition_).norm() > options.positionTolerance_ ||
        std::abs(ga.receiver_->pathDistance_ - gb.receiver_->pathDistance_) > options.lengthTolerance_)
        return false;

    // Chord length avoids acos() and loss of precision for very small angles.
    const double directionTolerance = 2.0 * std::sin(0.5 * options.angleTolerance_);
    if ((a.launchDirection_ - b.launchDirection_).norm() > directionTolerance ||
        (ga.finalDirection_ - gb.finalDirection_).norm() > directionTolerance) return false;
    for (std::size_t i = 0; i < ga.reflections_.size(); ++i) {
        if (ga.reflections_[i].surfaceIndex_ != gb.reflections_[i].surfaceIndex_ ||
            (ga.reflections_[i].position_ - gb.reflections_[i].position_).norm() > options.positionTolerance_)
            return false;
    }
    return true;
}
}

RefinementResult refine_path(const Vec3 &transmitterPosition, const Vec3 &launchDirection,
    const Vec3 &receiverPosition, const std::vector<Surface> &surfaces,
    const std::vector<std::size_t> &surfaceSequence, const RefinementOptions &options) {
    if (options.maxIterations_ < 0 || options.maxBacktracks_ < 0 ||
        !positive_finite(options.receiverTolerance_) || !positive_finite(options.finiteDifferenceStep_) ||
        !positive_finite(options.maxDirectionStep_) || !positive_finite(options.cornerSeparationTolerance_) ||
        !std::isfinite(options.tMin_) || options.tMin_ < 0.0 ||
        std::isnan(options.maxPathDistance_) || options.maxPathDistance_ <= 0.0) {
        throw std::invalid_argument("Invalid path refinement options");
    }

    auto evaluate = [&](const Vec3 &direction) {
        return evaluate_sequence(transmitterPosition, direction, receiverPosition, surfaces, surfaceSequence, options.tMin_);
    };
    RefinementResult result;
    result.transmitterPosition_ = transmitterPosition;
    result.receiverPosition_ = receiverPosition;
    result.launchDirection_ = launchDirection;
    result.geometry_ = evaluate(launchDirection); // Also checks endpoint, direction, and sequence arguments.
    if (!usable(result.geometry_)) {
        result.status_ = result.geometry_.status_ == SequenceStatus::AmbiguousHit
            ? RefinementStatus::UnresolvedCorner : RefinementStatus::InvalidGeometry;
        return result;
    }

    bool ambiguousUpdate = false;
    auto failure = [&](RefinementStatus status) {
        result.status_ = ambiguousUpdate || collapsing_corner(result.geometry_, options.cornerSeparationTolerance_)
            ? RefinementStatus::UnresolvedCorner : status;
        return result;
    };

    while (true) {
        const double miss = result.geometry_.receiver_->missDistance_;
        if (miss <= options.receiverTolerance_) {
            // Re-evaluate the final direction before accepting the solution.
            result.geometry_ = evaluate(result.launchDirection_);
            if (!usable(result.geometry_) || result.geometry_.receiver_->missDistance_ > options.receiverTolerance_)
                result.status_ = RefinementStatus::InvalidGeometry;
            else if (result.geometry_.receiver_->pathDistance_ > options.maxPathDistance_)
                result.status_ = RefinementStatus::PathTooLong;
            else if (collapsing_corner(result.geometry_, options.cornerSeparationTolerance_))
                result.status_ = RefinementStatus::UnresolvedCorner;
            else
                result.status_ = RefinementStatus::Converged;
            return result;
        }
        if (result.iterations_ >= options.maxIterations_) {
            return failure(RefinementStatus::IterationLimit);
        }

        // Rebuild a local tangent frame at every accepted direction.
        const Vec3 u = result.launchDirection_.unitOrthogonal();
        const Vec3 v = result.launchDirection_.cross(u).normalized();
        const Vec3 residual = result.geometry_.receiver_->residual_;
        Eigen::Matrix<double, 3, 2> jacobian;
        for (int column = 0; column < 2; ++column) {
            const Vec3 axis = column == 0 ? u : v;
            double epsilon = options.finiteDifferenceStep_;
            bool found = false;
            for (int attempt = 0; attempt < 8; ++attempt) {
                const auto plus = evaluate((result.launchDirection_ + epsilon * axis).normalized());
                const auto minus = evaluate((result.launchDirection_ - epsilon * axis).normalized());
                if (usable(plus) && usable(minus)) {
                    jacobian.col(column) = (plus.receiver_->residual_ - minus.receiver_->residual_) / (2.0 * epsilon);
                    found = true;
                } else if (usable(plus)) {
                    jacobian.col(column) = (plus.receiver_->residual_ - residual) / epsilon;
                    found = true;
                } else if (usable(minus)) {
                    jacobian.col(column) = (residual - minus.receiver_->residual_) / epsilon;
                    found = true;
                }
                if (found) break;
                epsilon *= 0.5;
            }
            if (!found || !jacobian.col(column).allFinite()) {
                return failure(RefinementStatus::DerivativeUnavailable);
            }
        }

        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 3, 2>> qr(jacobian);
        qr.setThreshold(1e-10);
        if (qr.rank() < 2) {
            return failure(RefinementStatus::SingularJacobian);
        }
        Eigen::Vector2d step = qr.solve(-residual);
        if (!step.allFinite() || step.norm() == 0.0) {
            return failure(RefinementStatus::Stalled);
        }
        if (step.norm() > options.maxDirectionStep_)
            step *= options.maxDirectionStep_ / step.norm();

        bool accepted = false;
        double scale = 1.0;
        for (int backtrack = 0; backtrack <= options.maxBacktracks_; ++backtrack) {
            const Vec3 trialDirection = (result.launchDirection_ + scale * (step.x() * u + step.y() * v)).normalized();
            auto trial = evaluate(trialDirection);
            ambiguousUpdate = ambiguousUpdate || trial.status_ == SequenceStatus::AmbiguousHit;
            if (usable(trial) && trial.receiver_->missDistance_ < miss) {
                result.launchDirection_ = trialDirection;
                result.geometry_ = std::move(trial);
                ++result.iterations_;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if (!accepted) {
            return failure(RefinementStatus::Stalled);
        }
    }
}

const char *refinement_status_name(RefinementStatus status) {
    switch (status) {
        case RefinementStatus::Converged: return "converged";
        case RefinementStatus::InvalidGeometry: return "invalid geometry or visibility";
        case RefinementStatus::DerivativeUnavailable: return "derivative unavailable";
        case RefinementStatus::SingularJacobian: return "singular Jacobian";
        case RefinementStatus::Stalled: return "stalled";
        case RefinementStatus::IterationLimit: return "iteration limit";
        case RefinementStatus::PathTooLong: return "path exceeds distance limit";
        case RefinementStatus::UnresolvedCorner: return "unresolved corner (retain received candidate)";
    }
    return "unknown";
}

std::vector<std::size_t> deduplicate_paths(const std::vector<RefinementResult> &results,
    const DeduplicationOptions &options) {
    if (!positive_finite(options.positionTolerance_) || !positive_finite(options.angleTolerance_) ||
        options.angleTolerance_ > std::numbers::pi || !positive_finite(options.lengthTolerance_)) {
        throw std::invalid_argument("Invalid path deduplication tolerances");
    }
    std::vector<std::size_t> candidates;
    for (std::size_t i = 0; i < results.size(); ++i)
        if (results[i].status_ == RefinementStatus::Converged && usable(results[i].geometry_))
            candidates.push_back(i);
    std::stable_sort(candidates.begin(), candidates.end(), [&](std::size_t a, std::size_t b) {
        return results[a].geometry_.receiver_->missDistance_ < results[b].geometry_.receiver_->missDistance_;
    });

    std::vector<std::size_t> unique;
    for (std::size_t candidate : candidates) {
        const bool duplicate = std::any_of(unique.begin(), unique.end(), [&](std::size_t representative) {
            return same_path(results[candidate], results[representative], options);
        });
        if (!duplicate) unique.push_back(candidate);
    }
    return unique;
}
