#include "Spreading.hpp"

#include <Eigen/SVD>
#include <cmath>
#include <stdexcept>

SpreadingResult calculate_spreading(
    const Vec3 &transmitterPosition, const Vec3 &launchDirection,
    const Vec3 &receiverPosition, const std::vector<Surface> &surfaces,
    const std::vector<std::size_t> &surfaceSequence, const SpreadingOptions &options) {
    const auto positive = [](double x) { return std::isfinite(x) && x > 0.0; };
    if (!positive(options.angularStep_) || options.angularStep_ > 0.1 ||
        options.maxStepHalvings_ < 0 || !positive(options.receiverTolerance_) ||
        !positive(options.relativeDerivativeTolerance_) ||
        !positive(options.absoluteDerivativeTolerance_) ||
        !positive(options.minimumSingularValue_) ||
        !positive(options.minimumSingularValueRatio_) || options.minimumSingularValueRatio_ >= 1.0 ||
        !positive(options.referenceDistance_) || !std::isfinite(options.tMin_) || options.tMin_ < 0.0)
        throw std::invalid_argument("Invalid spreading options");

    const auto central = evaluate_sequence(transmitterPosition, launchDirection,
                                           receiverPosition, surfaces, surfaceSequence, options.tMin_);
    SpreadingResult result;
    if (central.status_ != SequenceStatus::Valid || !central.receiver_ ||
        central.receiver_->missDistance_ > options.receiverTolerance_)
        return result;

    result.launchU_ = launchDirection.unitOrthogonal();
    result.launchV_ = launchDirection.cross(result.launchU_).normalized();
    const Vec3 planeNormal = central.finalDirection_;
    result.receiverU_ = planeNormal.unitOrthogonal();
    result.receiverV_ = planeNormal.cross(result.receiverU_).normalized();

    const auto coordinates = [&](double alpha, double beta) -> std::optional<Eigen::Vector2d> {
        const Vec3 direction = (launchDirection + alpha * result.launchU_ + beta * result.launchV_).normalized();
        const auto trial = evaluate_sequence(transmitterPosition, direction, receiverPosition,
                                             surfaces, surfaceSequence, options.tMin_);
        // The closest-point visibility interval may differ from the fixed-plane
        // interval. Recheck visibility below against the actual plane distance.
        if ((trial.status_ != SequenceStatus::Valid && trial.status_ != SequenceStatus::FinalSegmentBlocked) ||
            trial.reflections_.size() != surfaceSequence.size())
            return std::nullopt;
        const double denominator = trial.finalDirection_.dot(planeNormal);
        if (denominator <= 1e-8) return std::nullopt;
        const double distance = (receiverPosition - trial.finalOrigin_).dot(planeNormal) / denominator;
        if (!std::isfinite(distance) || distance < 0.0) return std::nullopt;
        if (distance > options.tMin_ && nearest_surface(trial.finalOrigin_, trial.finalDirection_,
                                                        surfaces, options.tMin_, distance))
            return std::nullopt;
        const Vec3 offset = trial.finalOrigin_ + distance * trial.finalDirection_ - receiverPosition;
        Eigen::Vector2d value{result.receiverU_.dot(offset), result.receiverV_.dot(offset)};
        if (!value.allFinite()) return std::nullopt;
        return value;
    };
    const auto derivative = [&](double h) -> std::optional<Eigen::Matrix2d> {
        const auto pu = coordinates(h, 0), mu = coordinates(-h, 0);
        const auto pv = coordinates(0, h), mv = coordinates(0, -h);
        if (!pu || !mu || !pv || !mv) return std::nullopt;
        Eigen::Matrix2d j;
        j.col(0) = (*pu - *mu) / (2.0 * h);
        j.col(1) = (*pv - *mv) / (2.0 * h);
        if (!j.allFinite()) return std::nullopt;
        return j;
    };

    double h = options.angularStep_;
    result.status_ = SpreadingStatus::InvalidPerturbation;
    for (int attempt = 0; attempt <= options.maxStepHalvings_; ++attempt, h *= 0.5) {
        const auto coarse = derivative(h), fine = derivative(h * 0.5);
        if (!coarse || !fine) continue;
        result.status_ = SpreadingStatus::UnstableDerivative;
        result.jacobian_ = *fine;
        result.angularStep_ = h * 0.5;
        result.derivativeDifference_ = (*coarse - *fine).norm();
        result.singularValues_ = Eigen::JacobiSVD<Eigen::Matrix2d>(*fine).singularValues();
        result.areaPerSolidAngle_ = std::abs(fine->determinant());
        if (result.derivativeDifference_ > options.absoluteDerivativeTolerance_ +
            options.relativeDerivativeTolerance_ * fine->norm())
            continue;
        // A large singular value must not hide uncertainty in a small one.
        const Eigen::Vector2d coarseSingularValues = Eigen::JacobiSVD<Eigen::Matrix2d>(*coarse).singularValues();
        bool stable = true;
        for (int i = 0; i < 2; ++i)
            stable = stable && std::abs(coarseSingularValues[i] - result.singularValues_[i]) <=
                     options.absoluteDerivativeTolerance_ + options.relativeDerivativeTolerance_ * result.
                     singularValues_[i];
        if (!stable) continue;
        if (result.singularValues_[1] <= options.minimumSingularValue_ ||
            result.singularValues_[1] <= options.minimumSingularValueRatio_ * result.singularValues_[0]) {
            result.status_ = SpreadingStatus::Caustic;
            return result;
        }
        const double factor = options.referenceDistance_ / std::sqrt(result.areaPerSolidAngle_);
        if (!std::isfinite(result.areaPerSolidAngle_) || !positive(factor)) return result;
        result.fieldFactor_ = factor;
        result.status_ = SpreadingStatus::Valid;
        return result;
    }
    return result;
}

const char *spreading_status_name(SpreadingStatus status) {
    switch (status) {
        case SpreadingStatus::Valid: return "valid";
        case SpreadingStatus::InvalidCentralPath: return "invalid central path";
        case SpreadingStatus::InvalidPerturbation: return "invalid neighbouring ray";
        case SpreadingStatus::UnstableDerivative: return "unstable derivative";
        case SpreadingStatus::Caustic: return "near receiver caustic";
    }
    return "unknown";
}
