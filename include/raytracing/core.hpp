#pragma once
#include <Eigen/Dense>
#include <complex>
#include <numbers>

// SI units throughout; phasors use the exp(-i omega t) convention.
namespace rt {
    using Vec3 = Eigen::Vector3d;
    using Vec3C = Eigen::Vector3cd;
    using Complex = std::complex<double>;

    inline constexpr double pi = std::numbers::pi;
    inline constexpr double speedOfLight = 299792458.0; // m/s, exact
    inline constexpr double vacuumPermeability = 1.25663706212e-6; // H/m, CODATA 2018
    inline constexpr double vacuumPermittivity = 1.0 / (vacuumPermeability * speedOfLight * speedOfLight);
    inline constexpr double freeSpaceImpedance = vacuumPermeability * speedOfLight; // ohm

    inline double wavenumber(double frequencyHz) { return 2.0 * pi * frequencyHz / speedOfLight; }
} // namespace rt
