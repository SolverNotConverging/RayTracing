#pragma once
#include "geometry.hpp"
#include <functional>

// Coarse geometric ray tracing used to discover candidate reflection sequences.
// Rays carry no field: fields are rebuilt from refined point-to-point paths.
namespace rt {
    struct ReflectionEvent {
        std::size_t surfaceIndex;
        std::size_t vertexIndex; // index of the hit point in TracedRay::vertices
        Vec3 normal;
    };

    struct TracedRay {
        std::vector<Vec3> vertices; // launch point, then every hit/stop point
        std::vector<ReflectionEvent> reflections;
        Vec3 direction; // current unit direction
        double distance = 0.0; // physical length travelled

        TracedRay(const Vec3 &origin, const Vec3 &launchDirection);
        const Vec3 &position() const { return vertices.back(); }
        void advance(double step);
        void reflect(const Vec3 &unitNormal, std::size_t surfaceIndex);
    };

    Vec3 reflected_direction(const Vec3 &incidentDirection, const Vec3 &unitNormal);

    // Entry distance into the reception sphere, or 0 when already inside.
    std::optional<double> intersect_receiver(const TracedRay &ray, const Vec3 &receiverPosition, double radius);

    enum class TraceStatus { Received, Escaped, ReflectionLimit, DistanceLimit, AmbiguousHit };
    inline constexpr std::size_t traceStatusCount = 5;

    struct TraceOptions {
        int maxReflections = 8;
        double maxDistance = 20.0;
    };

    struct TraceResult {
        TraceStatus status;
        int reflections;
        double distance;
    };

    // Trace one ray. Without a callback the reception sphere absorbs the ray.
    // With onReception, every segment that enters the sphere is reported (as a
    // snapshot ending at the sphere entry) and tracing continues to the physical
    // termination. Colocated Tx/Rx endpoints only report reflected returns.
    TraceResult trace_ray(TracedRay &ray, const std::vector<Surface> &surfaces, const Vec3 &receiverPosition,
                          double receiverRadius, const TraceOptions &options = {},
                          const std::function<void(const TracedRay &)> &onReception = {});

    // Quasi-uniform spherical Fibonacci launch directions.
    std::vector<Vec3> launch_directions(std::size_t count);
} // namespace rt
