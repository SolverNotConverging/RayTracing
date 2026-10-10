#pragma once
#include <numbers>

// Physical constants in SI units. The speed of light is exact by definition of
// the metre; mu0 is the CODATA 2018 recommended value, and eps0/eta0 follow
// from c^2 = 1/(mu0 eps0) and eta0 = mu0 c.
namespace rt::constants {
    inline constexpr double pi = std::numbers::pi;
    inline constexpr double speedOfLight = 299792458.0; // m/s
    inline constexpr double vacuumPermeability = 1.25663706212e-6; // H/m
    inline constexpr double vacuumPermittivity =
            1.0 / (vacuumPermeability * speedOfLight * speedOfLight); // F/m
    inline constexpr double freeSpaceImpedance = vacuumPermeability * speedOfLight; // ohm
} // namespace rt::constants
