#include <raytracing/sequence.hpp>

#include <algorithm>

namespace rt {
SequenceEvaluation evaluate_sequence(
    const Vec3 &transmitterPosition,
    const Vec3 &launchDirection,
    const Vec3 &receiverPosition,
    const std::vector<Surface> &surfaces,
    const std::vector<std::size_t> &surfaceSequence,
    double tMin) {
    SequenceEvaluation result;
    result.finalOrigin = transmitterPosition;
    result.finalDirection = launchDirection;
    result.reflections.reserve(surfaceSequence.size());
    double reflectedDistance = 0.0;

    for (std::size_t expectedIndex: surfaceSequence) {
        // Query the full scene: intersecting only the expected object could
        // incorrectly let a trial pass through another reflector.
        const auto hit = nearest_surface(result.finalOrigin, result.finalDirection, surfaces, tMin);
        if (!hit) {
            result.status = SequenceStatus::MissedSurface;
            return result;
        }
        if (hit->ambiguous) {
            result.finalOrigin += hit->distance * result.finalDirection;
            result.status = SequenceStatus::AmbiguousHit;
            result.blockingSurface = hit->surfaceIndex;
            return result;
        }
        // Equivalent coplanar faces can share an edge; accept either identity.
        if (std::find(hit->coincidentSurfaces.begin(), hit->coincidentSurfaces.end(), expectedIndex) ==
            hit->coincidentSurfaces.end()) {
            result.status = SequenceStatus::UnexpectedSurface;
            result.blockingSurface = hit->surfaceIndex;
            return result;
        }

        result.finalOrigin += hit->distance * result.finalDirection;
        reflectedDistance += hit->distance;
        result.reflections.push_back({expectedIndex, result.finalOrigin, hit->normal, hit->distance});
        result.finalDirection = (result.finalDirection -
                                  2.0 * result.finalDirection.dot(hit->normal) * hit->normal).normalized();
    }

    // Orthogonal projection onto the final ray's supporting line.
    const double finalDistance = result.finalDirection.dot(receiverPosition - result.finalOrigin);
    if (finalDistance < 0.0) {
        result.status = SequenceStatus::ReceiverBehind;
        return result;
    }
    const Vec3 closestPoint = result.finalOrigin + finalDistance * result.finalDirection;
    const Vec3 residual = closestPoint - receiverPosition;
    result.receiver = ReceiverApproach{
        closestPoint, residual, residual.norm(), finalDistance, reflectedDistance + finalDistance
    };

    // Visibility is checked along the actual trial ray, up to its closest point.
    // An off-axis trial is not connected to Rx by an artificial extra segment.
    if (finalDistance > tMin) {
        const auto blocker = nearest_surface(result.finalOrigin, result.finalDirection,
                                             surfaces, tMin, finalDistance);
        if (blocker) {
            result.status = blocker->ambiguous ? SequenceStatus::AmbiguousHit : SequenceStatus::FinalSegmentBlocked;
            result.blockingSurface = blocker->surfaceIndex;
        }
    }
    return result;
}
} // namespace rt
