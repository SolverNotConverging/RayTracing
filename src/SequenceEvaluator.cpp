#include "SequenceEvaluator.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>

SequenceEvaluation evaluate_sequence(
    const Vec3 &transmitterPosition,
    const Vec3 &launchDirection,
    const Vec3 &receiverPosition,
    const std::vector<Surface> &surfaces,
    const std::vector<std::size_t> &surfaceSequence,
    double tMin) {
    if (!transmitterPosition.allFinite() || !receiverPosition.allFinite() ||
        !launchDirection.allFinite() || std::abs(launchDirection.norm() - 1.0) > 1e-10 ||
        !std::isfinite(tMin) || tMin < 0.0) {
        throw std::invalid_argument("Sequence evaluation requires finite endpoints, a unit direction, and finite nonnegative tMin");
    }
    for (std::size_t index : surfaceSequence) {
        if (index >= surfaces.size()) {
            throw std::invalid_argument("Reflection sequence contains an invalid surface index");
        }
    }

    SequenceEvaluation result;
    result.finalOrigin_ = transmitterPosition;
    result.finalDirection_ = launchDirection;
    result.reflections_.reserve(surfaceSequence.size());
    double reflectedDistance = 0.0;

    for (std::size_t expectedIndex : surfaceSequence) {
        // Query the full scene: intersecting only the expected object could
        // incorrectly let a trial pass through another reflector.
        const auto hit = nearest_surface(result.finalOrigin_, result.finalDirection_, surfaces, tMin);
        if (!hit) {
            result.status_ = SequenceStatus::MissedSurface;
            return result;
        }
        if (hit->ambiguous_) {
            result.finalOrigin_ += hit->distance_ * result.finalDirection_;
            result.status_ = SequenceStatus::AmbiguousHit;
            result.blockingSurfaceIndex_ = hit->surfaceIndex_;
            return result;
        }
        // Equivalent coplanar faces can share an edge; accept either identity.
        if (std::find(hit->coincidentSurfaceIndices_.begin(), hit->coincidentSurfaceIndices_.end(), expectedIndex)
                == hit->coincidentSurfaceIndices_.end()) {
            result.status_ = SequenceStatus::UnexpectedSurface;
            result.blockingSurfaceIndex_ = hit->surfaceIndex_;
            return result;
        }

        result.finalOrigin_ += hit->distance_ * result.finalDirection_;
        reflectedDistance += hit->distance_;
        result.reflections_.push_back({expectedIndex, result.finalOrigin_, hit->normal_, hit->distance_});
        result.finalDirection_ = (result.finalDirection_ -
            2.0 * result.finalDirection_.dot(hit->normal_) * hit->normal_).normalized();
    }

    // Orthogonal projection onto the final ray's supporting line.
    const double finalDistance = result.finalDirection_.dot(receiverPosition - result.finalOrigin_);
    if (finalDistance < 0.0) {
        result.status_ = SequenceStatus::ReceiverBehind;
        return result;
    }
    const Vec3 closestPoint = result.finalOrigin_ + finalDistance * result.finalDirection_;
    const Vec3 residual = closestPoint - receiverPosition;
    result.receiver_ = ReceiverApproach{
        closestPoint, residual, residual.norm(), finalDistance, reflectedDistance + finalDistance
    };

    // Visibility is checked along the actual trial ray, up to its closest point.
    // An off-axis trial is not connected to Rx by an artificial extra segment.
    if (finalDistance > tMin) {
        const auto blocker = nearest_surface(result.finalOrigin_, result.finalDirection_,
                                              surfaces, tMin, finalDistance);
        if (blocker) {
            result.status_ = blocker->ambiguous_ ? SequenceStatus::AmbiguousHit : SequenceStatus::FinalSegmentBlocked;
            result.blockingSurfaceIndex_ = blocker->surfaceIndex_;
        }
    }
    return result;
}
