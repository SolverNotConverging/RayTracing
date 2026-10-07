#include "EMRay.hpp"

namespace {
void require_unit_normal(const Vec3 &normal) {
    if (!normal.allFinite() || std::abs(normal.norm() - 1.0) > 1e-10) {
        throw std::invalid_argument("Plane normal must be a unit vector");
    }
}
}

std::optional<double> intersect_plane(
    const EMRay &ray,
    const Vec3 &planePoint,
    const Vec3 &planeNormal) {
    require_unit_normal(planeNormal);
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


void EMRay::propagate(
    double distance,
    double refractiveIndex) {
    const double deltaOPL =
            refractiveIndex * distance;

    const double k0 =
            2.0 * std::numbers::pi / wavelength0_;

    const double phase =
            k0 * deltaOPL;

    const Complex phaseFactor =
            std::exp(Complex{0.0, phase});

    // Store the position, field phasor, and optical path for this step together.
    const RaySample &previous = path_.back();
    path_.push_back({
        previous.position + distance * direction_,
        previous.E * phaseFactor,
        previous.opticalPath + deltaOPL
    });
}

Vec3 reflected_direction(
    const Vec3 &incidentDirection,
    const Vec3 &unitNormal) {
    require_unit_normal(unitNormal);
    return incidentDirection
           - 2.0 * incidentDirection.dot(unitNormal) * unitNormal;
}

void EMRay::reflect_pec(const Vec3 &unitNormal) {
    const Vec3 reflected = reflected_direction(direction_, unitNormal);
    const RaySample &incident = path_.back();
    const Vec3C normal = unitNormal.cast<Complex>();

    // At a PEC, the reflected tangential E cancels the incident tangential E.
    const Vec3C reflectedE = -incident.E + 2.0 * normal.dot(incident.E) * normal;
    path_.push_back({incident.position, reflectedE, incident.opticalPath});
    direction_ = reflected;
}
