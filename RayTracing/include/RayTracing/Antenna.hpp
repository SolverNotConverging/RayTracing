#pragma once
#include "Material.hpp"
#include <filesystem>
#include <memory>
#include <string>

namespace rt {
    using ::Complex;
    using ::Cylinder;
    using ::Disk;
    using ::Rectangle;
    using ::Sphere;
    using ::Surface;
    using ::Triangle;
    using ::Vec3;
    using ::Vec3C;

    enum class Polarization { Vertical, Horizontal, RightCircular, LeftCircular };

    enum class AntennaKind {
        Isotropic,
        ShortDipole,
        ThinWireDipole,
        RectangularAperture,
        Imported
    };

    enum class PhasorConvention { PositiveTime, NegativeTime };

    enum class PowerReference { Accepted, Radiated, Stimulated };

    struct Isotropic {
        Polarization polarization = Polarization::Vertical;
        Complex amplitude = 1.0;
    };

    struct ShortDipole {
        double effectiveLengthMetres = 1e-4; // Uniform-current element, local z.
        Complex currentAmperes = 1.0;
    };

    struct ThinWireDipole {
        double lengthMetres = 0.001947; // Sinusoidal current, local z; centre-fed.
        Complex feedCurrentAmperes = 1.0;
    };

    struct RectangularAperture {
        double widthMetres = 0.01, heightMetres = 0.01;
        Complex fieldX = 1.0, fieldY = 0.0; // Uniform tangential E, V/m.
        bool pecBacked = true; // Infinite PEC baffle; local +z radiation hemisphere.
    };

    struct FarfieldData {
        std::vector<double> frequenciesHz, thetaDegrees, phiDegrees;
        // Cartesian rE vectors (volts), in antenna-local axes;
        // (f*phiCount+phi)*thetaCount+theta.
        std::vector<Vec3C> coefficients;
        std::vector<double> radiatedPower, acceptedPower, stimulatedPower;
        Vec3 exportPosition = Vec3::Zero();
        Eigen::Matrix3d exportOrientation = Eigen::Matrix3d::Identity();
        std::string provenance;
        bool frequencyIndependent = false;
        PhasorConvention inputConvention = PhasorConvention::PositiveTime;
        double coefficientScale = 1.0;
    };

    struct FarfieldImportOptions {
        PhasorConvention inputConvention = PhasorConvention::PositiveTime;
        PowerReference receivePowerReference = PowerReference::Accepted;
        double inputPowerWatts =
                0.0; // HFSS .ffd lacks this metadata; required for receive use.
        double coefficientScale = 1.0;
        bool honorExportAxes = true;
    };

    class Antenna {
    public:
        AntennaKind kind = AntennaKind::Isotropic;
        Eigen::Matrix3d orientation = Eigen::Matrix3d::Identity(); // Local to global.
        Isotropic isotropic;
        ShortDipole shortDipole;
        ThinWireDipole dipole;
        RectangularAperture aperture;
        std::shared_ptr<const FarfieldData> pattern;
        PowerReference receivePowerReference = PowerReference::Accepted;
        Complex receiveCalibration =
                1.0; // Explicit reciprocal receive port amplitude/phase calibration.
        Antenna() = default;

        Antenna(Isotropic model) : isotropic(model) {
        }

        Antenna(ShortDipole model)
            : kind(AntennaKind::ShortDipole), shortDipole(model) {
        }

        Antenna(ThinWireDipole model)
            : kind(AntennaKind::ThinWireDipole), dipole(model) {
        }

        Antenna(RectangularAperture model)
            : kind(AntennaKind::RectangularAperture), aperture(model) {
        }

        // Compose a right-hand rotation about a world-space axis (degrees).
        // Rotates both pattern directions and polarization; preserves export axes.
        Antenna &rotate(const Vec3 &worldAxis, double angleDegrees);

        void validate(double frequencyHz, const Medium &medium = {}) const;

        Vec3C farfield(const Vec3 &outwardDirection, double frequencyHz,
                       const Medium &medium = {}) const;

        // Normalized reciprocal port response in outward look direction. Contracts
        // bilinearly with incident E; time convention and Cartesian bases are shared.
        Complex receive(const Vec3 &lookDirection, const Vec3C &incidentField,
                        double frequencyHz, const Medium &medium = {}) const;

        double reference_power(double frequencyHz, const Medium &medium = {}) const;
    };

    Antenna load_farfield(const std::filesystem::path &path,
                          const FarfieldImportOptions &options = {});
} // namespace rt
