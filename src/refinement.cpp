#include <raytracing/refinement.hpp>

#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <utility>

namespace rt {
namespace {
    bool usable(const SequenceEvaluation &geometry) {
        return geometry.status == SequenceStatus::Valid && geometry.receiver &&
               geometry.receiver->residual.allFinite() &&
               std::isfinite(geometry.receiver->missDistance) &&
               std::isfinite(geometry.receiver->pathDistance);
    }

    bool collapsing_corner(const SequenceEvaluation &geometry, double tolerance) {
        for (std::size_t i = 1; i < geometry.reflections.size(); ++i) {
            const auto &a = geometry.reflections[i - 1];
            const auto &b = geometry.reflections[i];
            if ((a.position - b.position).norm() <= tolerance &&
                std::min((a.normal - b.normal).norm(), (a.normal + b.normal).norm()) > equivalentNormalTolerance)
                return true;
        }
        return false;
    }

    bool touches_cylinder_rim(const SequenceEvaluation &geometry, const std::vector<Surface> &surfaces,
                              double tolerance) {
        for (const auto &hit: geometry.reflections) {
            const auto *cylinder = std::get_if<Cylinder>(&surfaces.at(hit.surfaceIndex));
            if (!cylinder) continue;
            const Vec3 offset = hit.position - cylinder->center;
            const double axial = offset.dot(cylinder->axis);
            const double radial = (offset - axial * cylinder->axis).norm();
            // A rim has no smooth reflecting neighbourhood, including an open end.
            if (cylinder->halfLength - std::abs(axial) <= tolerance &&
                cylinder->radius - radial <= tolerance)
                return true;
        }
        return false;
    }

    bool same_path(const RefinementResult &a, const RefinementResult &b, const DeduplicationOptions &options) {
        const auto &ga = a.geometry;
        const auto &gb = b.geometry;
        if (ga.reflections.size() != gb.reflections.size()) return false;
        if ((a.transmitterPosition - b.transmitterPosition).norm() > options.positionTolerance ||
            (a.receiverPosition - b.receiverPosition).norm() > options.positionTolerance ||
            std::abs(ga.receiver->pathDistance - gb.receiver->pathDistance) > options.lengthTolerance)
            return false;

        // Chord length avoids acos() and loss of precision for very small angles.
        const double directionTolerance = 2.0 * std::sin(0.5 * options.angleTolerance);
        if ((a.launchDirection - b.launchDirection).norm() > directionTolerance ||
            (ga.finalDirection - gb.finalDirection).norm() > directionTolerance)
            return false;
        for (std::size_t i = 0; i < ga.reflections.size(); ++i) {
            if (ga.reflections[i].surfaceIndex != gb.reflections[i].surfaceIndex ||
                (ga.reflections[i].position - gb.reflections[i].position).norm() > options.positionTolerance)
                return false;
        }
        return true;
    }
}

RefinementResult refine_path(const Vec3 &transmitterPosition, const Vec3 &launchDirection,
                             const Vec3 &receiverPosition, const std::vector<Surface> &surfaces,
                             const std::vector<std::size_t> &surfaceSequence, const RefinementOptions &options) {
    auto evaluate = [&](const Vec3 &direction) {
        return evaluate_sequence(transmitterPosition, direction, receiverPosition, surfaces, surfaceSequence,
                                 options.tMin);
    };
    RefinementResult result;
    result.transmitterPosition = transmitterPosition;
    result.receiverPosition = receiverPosition;
    result.launchDirection = launchDirection;
    result.geometry = evaluate(launchDirection);
    if (!usable(result.geometry)) {
        result.status = result.geometry.status == SequenceStatus::AmbiguousHit
                             ? RefinementStatus::UnresolvedCorner
                             : RefinementStatus::InvalidGeometry;
        return result;
    }

    bool ambiguousUpdate = false;
    auto failure = [&](RefinementStatus status) {
        result.status = ambiguousUpdate || collapsing_corner(result.geometry, options.cornerSeparationTolerance)
                             ? RefinementStatus::UnresolvedCorner
                             : status;
        return result;
    };

    while (true) {
        const double miss = result.geometry.receiver->missDistance;
        if (miss <= options.receiverTolerance) {
            // Re-evaluate the final direction before accepting the solution.
            result.geometry = evaluate(result.launchDirection);
            if (!usable(result.geometry) || result.geometry.receiver->missDistance > options.receiverTolerance)
                result.status = RefinementStatus::InvalidGeometry;
            else if (result.geometry.receiver->pathDistance > options.maxPathDistance)
                result.status = RefinementStatus::PathTooLong;
            else if (collapsing_corner(result.geometry, options.cornerSeparationTolerance))
                result.status = RefinementStatus::UnresolvedCorner;
            else if (touches_cylinder_rim(result.geometry, surfaces, options.cornerSeparationTolerance))
                result.status = RefinementStatus::UnresolvedEdge;
            else
                result.status = RefinementStatus::Converged;
            return result;
        }
        if (result.iterations >= options.maxIterations) {
            return failure(RefinementStatus::IterationLimit);
        }

        // Rebuild a local tangent frame at every accepted direction.
        const Vec3 u = result.launchDirection.unitOrthogonal();
        const Vec3 v = result.launchDirection.cross(u).normalized();
        const Vec3 residual = result.geometry.receiver->residual;
        Eigen::Matrix<double, 3, 2> jacobian;
        for (int column = 0; column < 2; ++column) {
            const Vec3 axis = column == 0 ? u : v;
            double epsilon = options.finiteDifferenceStep;
            bool found = false;
            for (int attempt = 0; attempt < 8; ++attempt) {
                const auto plus = evaluate((result.launchDirection + epsilon * axis).normalized());
                const auto minus = evaluate((result.launchDirection - epsilon * axis).normalized());
                if (usable(plus) && usable(minus)) {
                    jacobian.col(column) = (plus.receiver->residual - minus.receiver->residual) / (2.0 * epsilon);
                    found = true;
                } else if (usable(plus)) {
                    jacobian.col(column) = (plus.receiver->residual - residual) / epsilon;
                    found = true;
                } else if (usable(minus)) {
                    jacobian.col(column) = (residual - minus.receiver->residual) / epsilon;
                    found = true;
                }
                if (found) break;
                epsilon *= 0.5;
            }
            if (!found || !jacobian.col(column).allFinite()) {
                return failure(RefinementStatus::DerivativeUnavailable);
            }
        }

        Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 3, 2> > qr(jacobian);
        qr.setThreshold(1e-10);
        if (qr.rank() < 2) {
            return failure(RefinementStatus::SingularJacobian);
        }
        Eigen::Vector2d step = qr.solve(-residual);
        if (!step.allFinite() || step.norm() == 0.0) {
            return failure(RefinementStatus::Stalled);
        }
        if (step.norm() > options.maxDirectionStep)
            step *= options.maxDirectionStep / step.norm();

        bool accepted = false;
        double scale = 1.0;
        for (int backtrack = 0; backtrack <= options.maxBacktracks; ++backtrack) {
            const Vec3 trialDirection = (result.launchDirection + scale * (step.x() * u + step.y() * v)).normalized();
            auto trial = evaluate(trialDirection);
            ambiguousUpdate = ambiguousUpdate || trial.status == SequenceStatus::AmbiguousHit;
            if (usable(trial) && trial.receiver->missDistance < miss) {
                result.launchDirection = trialDirection;
                result.geometry = std::move(trial);
                ++result.iterations;
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

std::vector<std::size_t> deduplicate_paths(const std::vector<RefinementResult> &results,
                                           const DeduplicationOptions &options) {
    std::vector<std::size_t> candidates;
    for (std::size_t i = 0; i < results.size(); ++i)
        if (results[i].status == RefinementStatus::Converged && usable(results[i].geometry))
            candidates.push_back(i);
    std::stable_sort(candidates.begin(), candidates.end(), [&](std::size_t a, std::size_t b) {
        return results[a].geometry.receiver->missDistance < results[b].geometry.receiver->missDistance;
    });

    std::vector<std::size_t> unique;
    for (std::size_t candidate: candidates) {
        const bool duplicate = std::any_of(unique.begin(), unique.end(), [&](std::size_t representative) {
            return same_path(results[candidate], results[representative], options);
        });
        if (!duplicate) unique.push_back(candidate);
    }
    return unique;
}
} // namespace rt
