#include "EMRay.hpp"

std::optional<double> intersect_plane(
    const EMRay &ray,
    const Vec3 &planePoint,
    const Vec3 &planeNormal) {
    const Vec3 &origin = ray.path_.back().position;

    const double denominator =
            planeNormal.dot(ray.direction_);

    constexpr double parallelTolerance = 1e-10;

    if (std::abs(denominator) < parallelTolerance) {
        return std::nullopt;
    }

    const double distance =
            planeNormal.dot(planePoint - origin) / denominator;

    constexpr double minimumDistance = 1e-8; // metres

    if (distance <= minimumDistance) {
        return std::nullopt;
    }

    return distance;
}


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
