#pragma once

#include "SequenceEvaluator.hpp"

struct RefinementOptions {
    int maxIterations_ = 30;
    double receiverTolerance_ = 1e-10; // metres; resolve focused roots before spreading/deduplication
    double finiteDifferenceStep_ = 1e-5; // tangent-direction perturbation, approximately radians
    double maxDirectionStep_ = 0.2; // norm of the tangent update
    int maxBacktracks_ = 12; // halve a rejected update this many times
    double tMin_ = 1e-8;
    double maxPathDistance_ = std::numeric_limits<double>::infinity();
    double cornerSeparationTolerance_ = 1e-5; // metres; unresolved corners and cylinder rims
};

enum class RefinementStatus {
    Converged,
    InvalidGeometry,
    DerivativeUnavailable,
    SingularJacobian,
    Stalled,
    IterationLimit,
    PathTooLong,
    UnresolvedCorner, // Preserve the coarse received candidate; not a verified point-Rx path.
    UnresolvedEdge // The stationary reflection lies on a cylinder rim.
};

struct RefinementResult {
    RefinementStatus status_ = RefinementStatus::IterationLimit;
    Vec3 transmitterPosition_;
    Vec3 receiverPosition_;
    Vec3 launchDirection_;
    SequenceEvaluation geometry_;
    int iterations_ = 0; // Number of accepted direction updates
};

// Local numerical shooting with a fixed ordered surface sequence.
// Geometry must be validated once before calling. No EM fields are rebuilt.
// Only Converged results passed the sequence/visibility, miss, and length checks.
RefinementResult refine_path(const Vec3 &transmitterPosition, const Vec3 &launchDirection,
                             const Vec3 &receiverPosition, const std::vector<Surface> &surfaces,
                             const std::vector<std::size_t> &surfaceSequence, const RefinementOptions &options = {});

const char *refinement_status_name(RefinementStatus status);

struct DeduplicationOptions {
    double positionTolerance_ = 1e-5; // metres, at Tx, Rx, and corresponding reflection points
    double angleTolerance_ = 1e-5; // radians, launch and arrival directions
    double lengthTolerance_ = 1e-5; // metres of physical path length
};

// Return indices into results, retaining only converged distinct paths.
// Representatives are considered in increasing receiver miss distance.
// Compare directly against representatives, never transitively merge clusters.
std::vector<std::size_t> deduplicate_paths(const std::vector<RefinementResult> &results,
                                           const DeduplicationOptions &options = {});
