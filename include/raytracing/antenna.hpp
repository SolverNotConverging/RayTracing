#pragma once
#include "material.hpp"
#include <filesystem>
#include <memory>
#include <string>
#include <variant>
#include <vector>

namespace rt {
    enum class Polarization { Vertical, Horizontal, RightCircular, LeftCircular };
    enum class PhasorConvention { PositiveTime, NegativeTime };
    enum class PowerReference { Accepted, Radiated, Stimulated };

    struct Isotropic {
        Polarization polarization = Polarization::Vertical;
        Complex amplitude = 1.0;
    };

    // Uniform-current element along local z.
    struct ShortDipole {
        double effectiveLength = 1e-4;
        Complex current = 1.0;
    };

    // Centre-fed sinusoidal-current wire along local z.
    struct ThinWireDipole {
        double length = 0.001947;
        Complex feedCurrent = 1.0;
    };

    // Uniform tangential E (V/m) on a PEC-backed aperture radiating into local +z.
    struct RectangularAperture {
        double width = 0.01, height = 0.01;
        Complex fieldX = 1.0, fieldY = 0.0;
    };

    // Sampled far field rE (volts) as Cartesian vectors in antenna-local axes,
    // indexed (frequency * phiCount + phi) * thetaCount + theta.
    struct FarfieldPattern {
        std::vector<double> frequenciesHz, thetaDegrees, phiDegrees;
        std::vector<Vec3C> coefficients;
        std::vector<double> radiatedPower, acceptedPower, stimulatedPower;
        Vec3 exportPosition = Vec3::Zero();
        Eigen::Matrix3d exportOrientation = Eigen::Matrix3d::Identity();
        bool frequencyIndependent = false;
        std::string source;
    };

    struct FarfieldImportOptions {
        PhasorConvention convention = PhasorConvention::PositiveTime;
        PowerReference receivePowerReference = PowerReference::Accepted;
        double inputPower = 0.0; // W; HFSS .ffd files lack it, but receive use requires it
        double scale = 1.0;
        bool useExportAxes = true;
    };

    using AntennaModel = std::variant<Isotropic, ShortDipole, ThinWireDipole, RectangularAperture,
                                      std::shared_ptr<const FarfieldPattern>>;

    class Antenna {
    public:
        AntennaModel model = Isotropic{};
        Eigen::Matrix3d orientation = Eigen::Matrix3d::Identity(); // local to global
        PowerReference receivePowerReference = PowerReference::Accepted;
        Complex receiveCalibration = 1.0;

        Antenna() = default;
        Antenna(AntennaModel antennaModel) : model(std::move(antennaModel)) {}
        Antenna(Isotropic m) : model(m) {}
        Antenna(ShortDipole m) : model(m) {}
        Antenna(ThinWireDipole m) : model(m) {}
        Antenna(RectangularAperture m) : model(m) {}

        // Right-hand rotation about a world axis, composed with the current orientation.
        Antenna &rotate(const Vec3 &worldAxis, double angleDegrees);

        // Throw std::invalid_argument for invalid parameters at this frequency.
        void validate(double frequencyHz, const Material &medium = {}) const;

        Vec3C farfield(const Vec3 &outwardDirection, double frequencyHz, const Material &medium = {}) const;

        // Reciprocal receive response for an incident field arriving from lookDirection
        // (outward from the antenna). Contracts bilinearly with the incident field.
        Complex receive(const Vec3 &lookDirection, const Vec3C &incidentField, double frequencyHz,
                        const Material &medium = {}) const;

        double reference_power(double frequencyHz, const Material &medium = {}) const;
    };

    // Load a CST .ffs (version 3.0) or HFSS .ffd far-field file.
    Antenna load_farfield(const std::filesystem::path &path, const FarfieldImportOptions &options = {});
} // namespace rt
