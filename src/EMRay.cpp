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
    if (!std::isfinite(distance) || distance < 0.0 ||
        !std::isfinite(refractiveIndex) || refractiveIndex <= 0.0) {
        throw std::invalid_argument(
            "Distance must be finite and nonnegative; refractive index must be finite and positive");
    }
    if (distance == 0.0) {
        return;
    }
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

std::optional<double> intersect_receiver(
    const EMRay &ray,
    const Vec3 &receiverPosition,
    double radius) {
    if (!receiverPosition.allFinite() || !std::isfinite(radius) || radius <= 0.0) {
        throw std::invalid_argument("Receiver position must be finite and radius must be finite and positive");
    }
    const Vec3 m =
            ray.path_.back().position - receiverPosition;

    if (m.norm() <= radius) {
        return 0.0;
    }

    // Project the sphere centre onto the ray's line.
    const double closestDistance = -ray.direction_.dot(m);
    if (closestDistance <= 0.0) {
        return std::nullopt;
    }

    // Compute the perpendicular separation directly instead of subtracting
    // two large squared distances in the quadratic discriminant.
    const Vec3 perpendicular = m + closestDistance * ray.direction_;
    const double halfChordSquared = radius * radius - perpendicular.squaredNorm();
    if (halfChordSquared < 0.0) {
        return std::nullopt;
    }

    return closestDistance - std::sqrt(halfChordSquared);
}

std::optional<double> intersect_rectangle(
    const EMRay &ray,
    const Rectangle &rectangle) {
    // Assumes orthonormal u and v, and positive half-sizes.
    const Vec3 normal = rectangle.u.cross(rectangle.v);

    const auto distance =
            intersect_plane(ray, rectangle.center, normal);

    if (!distance) {
        return std::nullopt;
    }

    const Vec3 hit =
            ray.path_.back().position + *distance * ray.direction_;

    const Vec3 offset = hit - rectangle.center;
    const double alpha = offset.dot(rectangle.u);
    const double beta = offset.dot(rectangle.v);

    constexpr double edgeTolerance = 1e-8; // metres

    if (std::abs(alpha) > rectangle.halfWidth + edgeTolerance ||
        std::abs(beta) > rectangle.halfHeight + edgeTolerance) {
        return std::nullopt;
    }

    return distance;
}

std::optional<SurfaceHit> nearest_surface(
    const EMRay &ray,
    const std::vector<Rectangle> &rectangles) {
    std::optional<SurfaceHit> nearest;

    for (std::size_t i = 0; i < rectangles.size(); ++i) {
        const auto distance =
                intersect_rectangle(ray, rectangles[i]);

        if (!distance) {
            continue;
        }

        if (!nearest || *distance < nearest->distance) {
            nearest = SurfaceHit{*distance, i};
        }
    }

    return nearest;
}

TraceResult trace_ray(EMRay &ray, const std::vector<Rectangle> &rectangles,
                      const Vec3 &receiverPosition, double receiverRadius,
                      const TraceOptions &options) {
    if (options.maxReflections < 0 || !std::isfinite(options.maxDistance) ||
        options.maxDistance <= 0.0 || !std::isfinite(options.refractiveIndex) ||
        options.refractiveIndex <= 0.0) {
        throw std::invalid_argument(
            "Trace limits and refractive index must be valid and positive (zero reflections is allowed)");
    }

    int reflections = 0;
    double travelled = 0.0;
    while (true) {
        const auto receiverDistance = intersect_receiver(ray, receiverPosition, receiverRadius);
        const auto surface = nearest_surface(ray, rectangles);
        const double remaining = options.maxDistance - travelled;
        const bool receiverFirst = receiverDistance &&
                                   (!surface || *receiverDistance < surface->distance);

        if (!receiverFirst && !surface) {
            // No future event: show a finite outgoing segment rather than an
            // infinite ray. This continuation still accumulates physical OPL.
            ray.propagate(remaining, options.refractiveIndex);
            return {TraceStatus::Escaped, reflections, options.maxDistance};
        }

        const double nextDistance = receiverFirst ? *receiverDistance : surface->distance;
        if (nextDistance > remaining) {
            ray.propagate(remaining, options.refractiveIndex);
            return {TraceStatus::DistanceLimit, reflections, options.maxDistance};
        }

        ray.propagate(nextDistance, options.refractiveIndex);
        travelled += nextDistance;
        if (receiverFirst) {
            return {TraceStatus::Received, reflections, travelled};
        }
        if (travelled >= options.maxDistance) {
            return {TraceStatus::DistanceLimit, reflections, travelled};
        }
        if (reflections >= options.maxReflections) {
            return {TraceStatus::ReflectionLimit, reflections, travelled};
        }

        const Rectangle &wall = rectangles[surface->rectangleIndex];
        ray.reflect_pec(wall.u.cross(wall.v));
        ++reflections;
    }
}

std::vector<Vec3> launch_directions(std::size_t count) {
    std::vector<Vec3> directions;
    directions.reserve(count);

    const double goldenAngle = PI * (3.0 - std::sqrt(5.0));

    for (std::size_t i = 0; i < count; ++i) {
        const double z =
                1.0 - 2.0 * (static_cast<double>(i) + 0.5)
                / static_cast<double>(count);

        const double radius = std::sqrt(1.0 - z * z);
        const double phi = goldenAngle * static_cast<double>(i);

        directions.emplace_back(
            radius * std::cos(phi),
            radius * std::sin(phi),
            z);
    }

    return directions;
}

Vec3C launch_polarization(const Vec3 &unitDirection)
{
    // Avoid projecting an axis nearly parallel to the ray.
    const Vec3 reference =
        std::abs(unitDirection.z()) < 0.9
        ? Vec3{0.0, 0.0, 1.0}
    : Vec3{1.0, 0.0, 0.0};

    const Vec3 transverse =
        reference - reference.dot(unitDirection) * unitDirection;

    return transverse.normalized().cast<Complex>();
}