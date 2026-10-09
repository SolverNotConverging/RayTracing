#include "FieldReconstruction.hpp"
#include "EMRay.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

ReconstructedField reconstruct_field(const RefinementResult &path,
                                     const SpreadingResult &spreading, const Vec3C &sourceReferenceField,
                                     const FieldReconstructionOptions &options) {
    const auto positive = [](double x) { return std::isfinite(x) && x > 0.0; };
    if (!positive(options.frequencyHz_) || !positive(options.refractiveIndex_) ||
        !positive(options.receiverTolerance_) || !positive(options.sourceReferenceAmplitude_) || !sourceReferenceField.
        allFinite() ||
        !path.launchDirection_.allFinite() ||
        std::abs(path.launchDirection_.norm() - 1.0) > 1e-10 ||
        std::abs(path.launchDirection_.cast<Complex>().dot(sourceReferenceField)) >
        1e-10 * sourceReferenceField.norm())
        throw std::invalid_argument(
            "Field reconstruction requires finite positive options and a transverse source field");
    const auto &geometry = path.geometry_;
    if (path.status_ != RefinementStatus::Converged || geometry.status_ != SequenceStatus::Valid ||
        !geometry.receiver_ || !geometry.receiver_->residual_.allFinite() ||
        !std::isfinite(geometry.receiver_->missDistance_) ||
        geometry.receiver_->missDistance_ > options.receiverTolerance_ ||
        geometry.receiver_->missDistance_ < 0 ||
        !positive(geometry.receiver_->pathDistance_) ||
        !std::isfinite(geometry.receiver_->finalSegmentDistance_) || geometry.receiver_->finalSegmentDistance_ < 0 ||
        spreading.status_ != SpreadingStatus::Valid || !spreading.fieldFactor_ ||
        !positive(*spreading.fieldFactor_))
        throw std::invalid_argument("Field reconstruction requires a converged receiver path and valid spreading");

    Vec3C field = sourceReferenceField;
    Vec3 direction = path.launchDirection_;
    double length = geometry.receiver_->finalSegmentDistance_;
    for (const auto &hit: geometry.reflections_) {
        if (!hit.normal_.allFinite() || std::abs(hit.normal_.norm() - 1.0) > 1e-10 ||
            !std::isfinite(hit.segmentDistance_) || hit.segmentDistance_ <= 0)
            throw std::invalid_argument("Invalid reflection in refined geometry");
        const Vec3C normal = hit.normal_.cast<Complex>();
        field = (-field + 2.0 * normal.dot(field) * normal).eval();
        direction = reflected_direction(direction, hit.normal_).normalized();
        length += hit.segmentDistance_;
    }
    if (!geometry.finalDirection_.allFinite() || (direction - geometry.finalDirection_).norm() > 1e-8 ||
        !std::isfinite(length) || std::abs(length - geometry.receiver_->pathDistance_) >
        1e-10 * std::max(1.0, length))
        throw std::invalid_argument("Inconsistent refined path length or arrival direction");

    ReconstructedField result;
    result.pathDistance_ = length;
    result.opticalPath_ = options.refractiveIndex_ * length;
    result.delaySeconds_ = result.opticalPath_ / C;
    result.frequencyHz_ = options.frequencyHz_;
    result.fieldFactor_ = *spreading.fieldFactor_;
    result.arrivalDirection_ = direction;
    result.transportedReferenceField_ = field;
    const double cycles = options.frequencyHz_ * result.delaySeconds_;
    if (!std::isfinite(result.opticalPath_) || !std::isfinite(cycles))
        throw std::invalid_argument("Optical path or propagation phase overflow");
    const Complex phase = std::polar(1.0, 2.0 * PI * std::remainder(cycles, 1.0));
    result.receiverField_ = result.fieldFactor_ * phase * field;
    result.normalizedReceiverField_ = result.receiverField_ / options.sourceReferenceAmplitude_;
    if (!result.receiverField_.allFinite() || !result.normalizedReceiverField_.allFinite())
        throw std::invalid_argument("Reconstructed field overflow");
    return result;
}
