#include "Tracing.hpp"
#include "Constants.hpp"
#include <cmath>
#include <stdexcept>

namespace {
    void require_unit_normal(const Vec3 &normal) {
        if (!normal.allFinite() || std::abs(normal.norm() - 1.0) > 1e-10)
            throw std::invalid_argument("Plane normal must be a unit vector");
    }
} // namespace

TracedRay::TracedRay(const Vec3 &position, const Vec3 &direction)
    : vertices_{position}, direction_(direction) {
    if (!position.allFinite())
        throw std::invalid_argument("Ray origin must be finite");
    if (!direction.allFinite() || direction.norm() == 0.0)
        throw std::invalid_argument("Direction must be finite and nonzero");
    direction_.normalize();
}

void TracedRay::advance(double distance) {
    if (!std::isfinite(distance) || distance < 0.0)
        throw std::invalid_argument("Distance must be finite and nonnegative");
    if (distance == 0.0)
        return;
    vertices_.push_back(position() + distance * direction_);
    distance_ += distance;
}

void TracedRay::reflect(const Vec3 &unitNormal, std::size_t surfaceIndex) {
    direction_ = reflected_direction(direction_, unitNormal);
    reflections_.push_back({surfaceIndex, vertices_.size() - 1, unitNormal});
}

Vec3 reflected_direction(const Vec3 &incidentDirection, const Vec3 &unitNormal) {
    require_unit_normal(unitNormal);
    return incidentDirection - 2.0 * incidentDirection.dot(unitNormal) * unitNormal;
}

std::optional<double> intersect_receiver(const TracedRay &ray, const Vec3 &receiverPosition,
                                         double radius) {
    if (!receiverPosition.allFinite() || !std::isfinite(radius) || radius <= 0.0)
        throw std::invalid_argument(
            "Receiver position must be finite and radius must be finite and positive");
    const Vec3 m = ray.position() - receiverPosition;
    if (m.norm() <= radius)
        return 0.0;

    // Project the sphere centre onto the ray's line.
    const double closestDistance = -ray.direction_.dot(m);
    if (closestDistance <= 0.0)
        return std::nullopt;

    // Compute the perpendicular separation directly instead of subtracting
    // two large squared distances in the quadratic discriminant.
    const Vec3 perpendicular = m + closestDistance * ray.direction_;
    const double halfChordSquared = radius * radius - perpendicular.squaredNorm();
    if (halfChordSquared < 0.0)
        return std::nullopt;
    return closestDistance - std::sqrt(halfChordSquared);
}

TraceResult trace_ray(TracedRay &ray, const std::vector<Surface> &surfaces,
                      const Vec3 &receiverPosition, double receiverRadius,
                      const TraceOptions &options,
                      const std::function<void(const TracedRay &)> &onReception) {
    if (options.maxReflections_ < 0 || !std::isfinite(options.maxDistance_) ||
        options.maxDistance_ <= 0.0)
        throw std::invalid_argument(
            "Trace limits must be valid and positive (zero reflections is allowed)");

    int reflections = 0;
    double travelled = 0.0;
    // A transmitter inside the reception sphere must first leave it. Otherwise
    // a monostatic trace terminates at launch instead of finding a return.
    bool receiverArmed = (ray.position() - receiverPosition).norm() > receiverRadius;
    // Validate the receiver even when reception starts disarmed.
    intersect_receiver(ray, receiverPosition, receiverRadius);
    while (true) {
        const auto receiverDistance = (receiverArmed || (onReception && reflections > 0))
                                          ? intersect_receiver(ray, receiverPosition, receiverRadius)
                                          : std::optional<double>{};
        const auto surface = nearest_surface(ray.position(), ray.direction_, surfaces);
        const double remaining = options.maxDistance_ - travelled;
        const bool receiverFirst = receiverDistance &&
                                   (!surface || *receiverDistance < surface->distance_);
        if (receiverFirst && onReception && *receiverDistance <= remaining) {
            TracedRay received = ray;
            received.advance(*receiverDistance);
            onReception(received);
        }
        const bool stopAtReceiver = receiverFirst && !onReception;

        if (!stopAtReceiver && !surface) {
            // No future event: show a finite outgoing segment rather than an
            // infinite ray.
            ray.advance(remaining);
            return {TraceStatus::Escaped, reflections, options.maxDistance_};
        }

        const double nextDistance = stopAtReceiver ? *receiverDistance : surface->distance_;
        if (nextDistance > remaining) {
            ray.advance(remaining);
            return {TraceStatus::DistanceLimit, reflections, options.maxDistance_};
        }

        ray.advance(nextDistance);
        travelled += nextDistance;
        if ((ray.position() - receiverPosition).norm() > receiverRadius)
            receiverArmed = true;
        if (stopAtReceiver)
            return {TraceStatus::Received, reflections, travelled};
        if (surface->ambiguous_)
            return {TraceStatus::AmbiguousHit, reflections, travelled};
        if (travelled >= options.maxDistance_)
            return {TraceStatus::DistanceLimit, reflections, travelled};
        if (reflections >= options.maxReflections_)
            return {TraceStatus::ReflectionLimit, reflections, travelled};

        ray.reflect(surface->normal_, surface->surfaceIndex_);
        ++reflections;
    }
}

std::vector<Vec3> launch_directions(std::size_t count) {
    std::vector<Vec3> directions;
    directions.reserve(count);
    const double goldenAngle = rt::constants::pi * (3.0 - std::sqrt(5.0));
    for (std::size_t i = 0; i < count; ++i) {
        const double z = 1.0 - 2.0 * (static_cast<double>(i) + 0.5) / static_cast<double>(count);
        const double radius = std::sqrt(1.0 - z * z);
        const double phi = goldenAngle * static_cast<double>(i);
        directions.emplace_back(radius * std::cos(phi), radius * std::sin(phi), z);
    }
    return directions;
}
