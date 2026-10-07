#include "EMRay.hpp"

void propagate(
    EMRay &ray,
    double distance,
    double refractiveIndex) {
    const double deltaOPL =
            refractiveIndex * distance;

    const double k0 =
            2.0 * std::numbers::pi / ray.wavelength0_;

    const double phase =
            k0 * deltaOPL;

    const Complex phaseFactor =
            std::exp(Complex{0.0, phase});

    // Store the position, field phasor, and optical path for this step together.
    const RaySample &previous = ray.path_.back();
    ray.path_.push_back({
        previous.position + distance * ray.direction_,
        previous.E * phaseFactor,
        previous.opticalPath + deltaOPL
    });
}
