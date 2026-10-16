#pragma once
#include "geometry.hpp"

namespace rt {
    enum class SequenceStatus {
        Valid, // correct reflection sequence and unobstructed approach; may still miss Rx
        MissedSurface,
        UnexpectedSurface,
        ReceiverBehind,
        FinalSegmentBlocked,
        AmbiguousHit
    };

    struct SequenceReflection {
        std::size_t surfaceIndex;
        Vec3 position;
        Vec3 normal;
        double segmentDistance; // from the previous vertex
    };

    struct ReceiverApproach {
        Vec3 closestPoint;
        Vec3 residual; // closestPoint - receiverPosition
        double missDistance;
        double finalSegmentDistance;
        double pathDistance; // total physical length to closestPoint
    };

    struct SequenceEvaluation {
        SequenceStatus status = SequenceStatus::Valid;
        std::vector<SequenceReflection> reflections; // replayed prefix on failure
        Vec3 finalOrigin = Vec3::Zero();
        Vec3 finalDirection = Vec3::UnitX();
        std::optional<ReceiverApproach> receiver; // after all reflections and a forward projection
        std::optional<std::size_t> blockingSurface;
    };

    // Replay exactly surfaceSequence with nearest hits and specular reflection.
    // An empty sequence evaluates the direct ray. Valid does not mean the
    // receiver was reached: check receiver->missDistance.
    SequenceEvaluation evaluate_sequence(const Vec3 &transmitterPosition, const Vec3 &launchDirection,
                                         const Vec3 &receiverPosition, const std::vector<Surface> &surfaces,
                                         const std::vector<std::size_t> &surfaceSequence, double tMin = 1e-8);
} // namespace rt
