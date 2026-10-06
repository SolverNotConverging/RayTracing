#pragma once

#include "crt/ray.hpp"

#include <cmath>
#include <complex>
#include <numbers>

namespace crt {

    constexpr double speed_of_light = 299'792'458.0;

    struct PathContribution {
        double length;
        double delay;

        std::complex<double> coefficient;
    };

    inline PathContribution freeSpacePath(
        const Vec3& tx,
        const Vec3& rx,
        double frequency)
    {
        const double length = (rx - tx).norm();

        const double delay =
            length / speed_of_light;

        const double wavelength =
            speed_of_light / frequency;

        const double k =
            2.0 * std::numbers::pi * frequency
            / speed_of_light;

        const double amplitude =
            wavelength /
            (4.0 * std::numbers::pi * length);

        const std::complex<double> phase =
            std::exp(
                std::complex<double>{
                    0.0,
                    -k * length
                }
            );

        return {
            .length = length,
            .delay = delay,
            .coefficient = amplitude * phase
        };
    }

} // namespace crt