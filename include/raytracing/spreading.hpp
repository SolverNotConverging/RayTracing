#pragma once
#include "sequence.hpp"

namespace rt {
    enum class SpreadingStatus { Valid, InvalidCentralPath, InvalidPerturbation, UnstableDerivative, Caustic };

    struct SpreadingOptions {
        double angularStep = 1e-4;
        int maxStepHalvings = 8;
        double receiverTolerance = 1e-10;
        double relativeDerivativeTolerance = 1e-5;
        double absoluteDerivativeTolerance = 1e-8; // metres per radian
        double minimumSingularValue = 1e-7; // metres; near-caustic cutoff
        double minimumSingularValueRatio = 1e-8;
        double referenceDistance = 1.0; // source field is specified on this sphere
        double tMin = 1e-8;
    };

    struct SpreadingResult {
        SpreadingStatus status = SpreadingStatus::InvalidCentralPath;
        // Rows: receiver-plane coordinates; columns: launch tangent coordinates.
        Eigen::Matrix2d jacobian = Eigen::Matrix2d::Zero();
        Eigen::Vector2d singularValues = Eigen::Vector2d::Zero(); // descending
        double areaPerSolidAngle = 0.0; // |det J|, m^2/sr
        double angularStep = 0.0;
        double derivativeDifference = 0.0; // ||J(h) - J(h/2)||_F
        std::optional<double> fieldFactor; // referenceDistance / sqrt(|det J|), Valid only
        unsigned int causticCount = 0; // internal conjugate points, with multiplicity
        Vec3 launchU = Vec3::Zero(), launchV = Vec3::Zero();
        Vec3 receiverU = Vec3::Zero(), receiverV = Vec3::Zero();
    };

    // Point-source geometrical spreading for a fixed specular sequence, by
    // central differences on a receiver plane, checked against half the step.
    SpreadingResult calculate_spreading(const Vec3 &transmitterPosition, const Vec3 &launchDirection,
                                        const Vec3 &receiverPosition, const std::vector<Surface> &surfaces,
                                        const std::vector<std::size_t> &surfaceSequence,
                                        const SpreadingOptions &options = {});
} // namespace rt
