#pragma once

#include "SequenceEvaluator.hpp"

enum class SpreadingStatus {
    Valid,
    InvalidCentralPath,
    InvalidPerturbation,
    UnstableDerivative,
    Caustic
};

struct SpreadingOptions {
    double angularStep_ = 1e-4;
    int maxStepHalvings_ = 8;
    double receiverTolerance_ = 1e-6; // Require a refined point-receiver path.
    double relativeDerivativeTolerance_ = 1e-5;
    double absoluteDerivativeTolerance_ = 1e-8; // metres per angular coordinate
    double minimumSingularValue_ = 1e-7; // metres; near-caustic cutoff
    double minimumSingularValueRatio_ = 1e-8;
    double referenceDistance_ = 1.0; // Field specified on this source-centred sphere.
    double tMin_ = 1e-8;
};

struct SpreadingResult {
    SpreadingStatus status_ = SpreadingStatus::InvalidCentralPath;
    // Rows: receiver tangent coordinates; columns: launch tangent coordinates.
    Eigen::Matrix2d jacobian_ = Eigen::Matrix2d::Zero();
    Eigen::Vector2d singularValues_ = Eigen::Vector2d::Zero(); // descending
    double areaPerSolidAngle_ = 0.0; // |det J|, m^2 / sr locally
    double angularStep_ = 0.0; // Step used for the returned derivative.
    double derivativeDifference_ = 0.0; // ||J(h) - J(h/2)||_F
    std::optional<double> fieldFactor_; // referenceDistance / sqrt(|det J|)
    Vec3 launchU_ = Vec3::Zero(), launchV_ = Vec3::Zero();
    Vec3 receiverU_ = Vec3::Zero(), receiverV_ = Vec3::Zero();
};

// Point-source geometrical spreading in a homogeneous medium, for a fixed
// specular reflection sequence. Validate surfaces first. Invalid arguments throw.
// Uses central differences on one fixed receiver plane and checks h against h/2.
// Only Valid supplies a field factor. No field transport or caustic phase is
// computed; a Valid endpoint does not exclude caustics earlier along the path.
SpreadingResult calculate_spreading(
    const Vec3 &transmitterPosition, const Vec3 &launchDirection,
    const Vec3 &receiverPosition, const std::vector<Surface> &surfaces,
    const std::vector<std::size_t> &surfaceSequence,
    const SpreadingOptions &options = {});

const char *spreading_status_name(SpreadingStatus status);
