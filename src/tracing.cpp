#include <raytracing/tracing.hpp>
#include <cmath>

namespace rt {
    TracedRay::TracedRay(const Vec3 &origin, const Vec3 &launchDirection)
        : vertices{origin}, direction(launchDirection.normalized()) {}

    void TracedRay::advance(double step) {
        if (step == 0.0) return;
        vertices.push_back(position() + step * direction);
        distance += step;
    }

    void TracedRay::reflect(const Vec3 &unitNormal, std::size_t surfaceIndex) {
        direction = reflected_direction(direction, unitNormal);
        reflections.push_back({surfaceIndex, vertices.size() - 1, unitNormal});
    }

    Vec3 reflected_direction(const Vec3 &incidentDirection, const Vec3 &unitNormal) {
        return incidentDirection - 2.0 * incidentDirection.dot(unitNormal) * unitNormal;
    }

    std::optional<double> intersect_receiver(const TracedRay &ray, const Vec3 &receiverPosition, double radius) {
        const Vec3 m = ray.position() - receiverPosition;
        if (m.norm() <= radius) return 0.0;
        const double closestDistance = -ray.direction.dot(m);
        if (closestDistance <= 0.0) return std::nullopt;
        // Perpendicular separation directly, avoiding cancellation in the discriminant.
        const Vec3 perpendicular = m + closestDistance * ray.direction;
        const double halfChordSquared = radius * radius - perpendicular.squaredNorm();
        if (halfChordSquared < 0.0) return std::nullopt;
        return closestDistance - std::sqrt(halfChordSquared);
    }

    TraceResult trace_ray(TracedRay &ray, const std::vector<Surface> &surfaces, const Vec3 &receiverPosition,
                          double receiverRadius, const TraceOptions &options,
                          const std::function<void(const TracedRay &)> &onReception) {
        int reflections = 0;
        double travelled = 0.0;
        // Colocated launches must leave the sphere before a direct reception counts.
        // Distinct endpoints may share a capture sphere; refinement validates them.
        const double endpointSeparation = (ray.position() - receiverPosition).norm();
        bool receiverArmed = endpointSeparation > receiverRadius || (onReception && endpointSeparation > 0.0);
        while (true) {
            const auto receiverDistance = receiverArmed || (onReception && reflections > 0)
                                              ? intersect_receiver(ray, receiverPosition, receiverRadius)
                                              : std::nullopt;
            const auto surface = nearest_surface(ray.position(), ray.direction, surfaces);
            const double remaining = options.maxDistance - travelled;
            const bool receiverFirst = receiverDistance && (!surface || *receiverDistance < surface->distance);
            if (receiverFirst && onReception && *receiverDistance <= remaining) {
                TracedRay received = ray;
                received.advance(*receiverDistance);
                onReception(received);
            }
            const bool stopAtReceiver = receiverFirst && !onReception;

            if (!stopAtReceiver && !surface) {
                ray.advance(remaining); // draw a finite outgoing segment
                return {TraceStatus::Escaped, reflections, options.maxDistance};
            }
            const double next = stopAtReceiver ? *receiverDistance : surface->distance;
            if (next > remaining) {
                ray.advance(remaining);
                return {TraceStatus::DistanceLimit, reflections, options.maxDistance};
            }

            ray.advance(next);
            travelled += next;
            if ((ray.position() - receiverPosition).norm() > receiverRadius) receiverArmed = true;
            if (stopAtReceiver) return {TraceStatus::Received, reflections, travelled};
            if (surface->ambiguous) return {TraceStatus::AmbiguousHit, reflections, travelled};
            if (travelled >= options.maxDistance) return {TraceStatus::DistanceLimit, reflections, travelled};
            if (reflections >= options.maxReflections) return {TraceStatus::ReflectionLimit, reflections, travelled};
            ray.reflect(surface->normal, surface->surfaceIndex);
            ++reflections;
        }
    }

    std::vector<Vec3> launch_directions(std::size_t count) {
        std::vector<Vec3> directions;
        directions.reserve(count);
        const double goldenAngle = pi * (3.0 - std::sqrt(5.0));
        for (std::size_t i = 0; i < count; ++i) {
            const double z = 1.0 - 2.0 * (static_cast<double>(i) + 0.5) / static_cast<double>(count);
            const double radius = std::sqrt(1.0 - z * z);
            const double phi = goldenAngle * static_cast<double>(i);
            directions.emplace_back(radius * std::cos(phi), radius * std::sin(phi), z);
        }
        return directions;
    }
} // namespace rt
