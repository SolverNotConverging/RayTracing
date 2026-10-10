#pragma once
#include "Surface.hpp"
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

// Coarse geometric ray tracing used to discover candidate reflection sequences.
// Rays carry no field: fields are rebuilt from refined point-to-point paths
// (FieldReconstruction), so the search is independent of antenna patterns and
// materials. Every surface reflects specularly here.

struct ReflectionEvent {
    std::size_t surfaceIndex_;
    std::size_t vertexIndex_; // Index of the hit point in TracedRay::vertices_
    Vec3 normal_;
};

struct TracedRay {
    std::vector<Vec3> vertices_; // Launch point, then every hit/stop point
    std::vector<ReflectionEvent> reflections_;
    Vec3 direction_; // Current unit direction
    double distance_ = 0.0; // Physical length travelled, metres

    // Throws unless position is finite and direction is finite and nonzero.
    TracedRay(const Vec3 &position, const Vec3 &direction);

    const Vec3 &position() const { return vertices_.back(); }

    void advance(double distance);

    // Specular reflection at the current position.
    void reflect(const Vec3 &unitNormal, std::size_t surfaceIndex);
};

Vec3 reflected_direction(const Vec3 &incidentDirection, const Vec3 &unitNormal);

// Returns the entry distance into the reception sphere, or 0 if already inside.
// receiverPosition must be finite and radius must be finite and positive.
std::optional<double> intersect_receiver(const TracedRay &ray, const Vec3 &receiverPosition,
                                         double radius);

enum class TraceStatus { Received, Escaped, ReflectionLimit, DistanceLimit, AmbiguousHit };

inline constexpr std::size_t traceStatusCount = 5;

struct TraceOptions {
    int maxReflections_ = 8;
    double maxDistance_ = 20.0; // physical distance in metres, from this call's start
};

struct TraceResult {
    TraceStatus status_;
    int reflections_;
    double distance_;
};

// Trace one ray. After the last permitted reflection, still check the next
// segment for reception. Escaping rays are drawn out to maxDistance_.
// When launched inside/on the Rx sphere, reception starts after leaving it.
// With onReception, report each segment intersecting the capture sphere and
// continue tracing to the physical termination. The sphere is then a candidate
// search region, not an absorbing object. Each callback receives a snapshot
// that ends at the sphere entry point.
TraceResult trace_ray(TracedRay &ray, const std::vector<Surface> &surfaces,
                      const Vec3 &receiverPosition, double receiverRadius,
                      const TraceOptions &options = {},
                      const std::function<void(const TracedRay &)> &onReception = {});

// Quasi-uniform spherical Fibonacci launch directions.
std::vector<Vec3> launch_directions(std::size_t count);
