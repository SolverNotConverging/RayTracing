#pragma once
#include "sequence.hpp"

namespace rt {
    struct RefinementOptions {
        int maxIterations = 30;
        double receiverTolerance = 1e-10; // metres
        double finiteDifferenceStep = 1e-5; // radians, approximately
        double maxDirectionStep = 0.2;
        int maxBacktracks = 12;
        double tMin = 1e-8;
        double maxPathDistance = infinity;
        double cornerSeparationTolerance = 1e-5; // metres
    };

    enum class RefinementStatus {
        Converged,
        InvalidGeometry,
        DerivativeUnavailable,
        SingularJacobian,
        Stalled,
        IterationLimit,
        PathTooLong,
        UnresolvedCorner, // coarse candidate kept; not a verified point-to-point path
        UnresolvedEdge // stationary reflection lies on a cylinder rim
    };

    struct RefinementResult {
        RefinementStatus status = RefinementStatus::IterationLimit;
        Vec3 transmitterPosition;
        Vec3 receiverPosition;
        Vec3 launchDirection;
        SequenceEvaluation geometry;
        int iterations = 0;
    };

    // Newton shooting on the launch direction for a fixed surface sequence.
    RefinementResult refine_path(const Vec3 &transmitterPosition, const Vec3 &launchDirection,
                                 const Vec3 &receiverPosition, const std::vector<Surface> &surfaces,
                                 const std::vector<std::size_t> &surfaceSequence,
                                 const RefinementOptions &options = {});

    struct DeduplicationOptions {
        double positionTolerance = 1e-5; // metres
        double angleTolerance = 1e-5; // radians
        double lengthTolerance = 1e-5; // metres
    };

    // Indices of distinct converged paths, best receiver miss first.
    std::vector<std::size_t> deduplicate_paths(const std::vector<RefinementResult> &paths,
                                               const DeduplicationOptions &options = {});
} // namespace rt
