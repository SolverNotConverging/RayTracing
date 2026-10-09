#pragma once

#include "PathRefinement.hpp"
#include "Spreading.hpp"

struct FieldReconstructionOptions {
    double frequencyHz_ = 77e9;
    double refractiveIndex_ = 1.0; // Homogeneous, nondispersive medium.
    double receiverTolerance_ = 1e-6;
    double sourceReferenceAmplitude_ = 1.0; // One common reference across all launch directions.
};

struct ReconstructedField {
    Vec3C receiverField_ = Vec3C::Zero();
    Vec3C normalizedReceiverField_ = Vec3C::Zero(); // Divide by the common excitation reference.
    Vec3C transportedReferenceField_ = Vec3C::Zero(); // PEC transport before phase/spreading.
    Vec3 arrivalDirection_ = Vec3::UnitX();
    double pathDistance_ = 0.0;
    double opticalPath_ = 0.0;
    double delaySeconds_ = 0.0;
    double frequencyHz_ = 0.0;
    double fieldFactor_ = 0.0;
};

// Rebuild the receiver field from the refined geometry, never from the coarse
// reception-sphere ray. The source field is transverse, specified at spreading's
// reference distance with a source phase reference at Tx. Uses exp(-i omega t),
// hence propagation exp(+i k0 opticalPath), consistent with EMRay::propagate.
// Requires Converged geometry and Valid spreading for that same path. Throws on
// invalid inputs. Caustic-crossing phase is not included.
ReconstructedField reconstruct_field(const RefinementResult &path,
                                     const SpreadingResult &spreading, const Vec3C &sourceReferenceField,
                                     const FieldReconstructionOptions &options = {});
