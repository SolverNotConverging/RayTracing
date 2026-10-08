#pragma once

#include "Surface.hpp"

enum class SequenceStatus {
    Valid, // Correct reflection sequence and unobstructed forward approach; may still miss Rx.
    MissedSurface,
    UnexpectedSurface,
    ReceiverBehind,
    FinalSegmentBlocked,
    AmbiguousHit
};

struct SequenceReflection {
    std::size_t surfaceIndex_;
    Vec3 position_;
    Vec3 normal_;
    double segmentDistance_; // Physical distance from the previous vertex
};

struct ReceiverApproach {
    Vec3 closestPoint_;
    Vec3 residual_; // closestPoint_ - receiverPosition
    double missDistance_;
    double finalSegmentDistance_; // Forward distance along the final ray
    double pathDistance_; // Total physical length to closestPoint_, not to sphere entry
};

struct SequenceEvaluation {
    SequenceStatus status_ = SequenceStatus::Valid;
    std::vector<SequenceReflection> reflections_; // Successfully replayed prefix on failure
    Vec3 finalOrigin_ = Vec3::Zero();
    Vec3 finalDirection_ = Vec3::UnitX();
    // Present only after all reflections and a forward receiver projection.
    // Also available on FinalSegmentBlocked for diagnostics.
    std::optional<ReceiverApproach> receiver_;
    // Populated for UnexpectedSurface, FinalSegmentBlocked, and AmbiguousHit.
    std::optional<std::size_t> blockingSurfaceIndex_;
};

// Replay exactly surfaceSequence using nearest hits and specular reflection.
// Surfaces must already be validated and remain unchanged during evaluation.
// Direction must be unit length. An empty sequence evaluates a direct ray.
// tMin excludes self-hits, in metres, just as in the surface intersection API.
// Invalid arguments throw; geometrically invalid trials return a failure status.
// Valid does NOT mean convergence: the caller must check receiver_->missDistance_.
// The receiver is a point; no reception sphere, field, wavelength, or phase is used.
SequenceEvaluation evaluate_sequence(
    const Vec3 &transmitterPosition,
    const Vec3 &launchDirection,
    const Vec3 &receiverPosition,
    const std::vector<Surface> &surfaces,
    const std::vector<std::size_t> &surfaceSequence,
    double tMin = 1e-8);
