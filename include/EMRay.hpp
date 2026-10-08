#pragma once
#include "Surface.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>
#include <optional>

constexpr double PI = std::numbers::pi;
constexpr double C = 3e8;

struct RaySample {
    Vec3 position_;
    Vec3C E_;
    double opticalPath_;
};

struct ReflectionEvent {
    std::size_t surfaceIndex_;
    std::size_t incidentSampleIndex_;
    Vec3 normal_;
};

struct EMRay {
    std::vector<RaySample> path_;
    std::vector<ReflectionEvent> reflections_;
    Vec3 direction_;

    double wavelength0_;

    EMRay(const Vec3 &position, const Vec3 &direction, const Vec3C &E, const double wavelength0) : path_{
            {position, E, 0.0}
        },
        direction_(direction),
        wavelength0_(wavelength0) {
        if (!direction.allFinite() || direction.norm() == 0.0) {
            throw std::invalid_argument("Direction must be finite and nonzero");
        }
        direction_.normalize();
        if (!std::isfinite(wavelength0) || wavelength0 <= 0.0) {
            throw std::invalid_argument("Wavelength must be finite and positive");
        }
        if (!E.allFinite() || E.norm() == 0.0) {
            throw std::invalid_argument("E must be finite and nonzero");
        }
        if (const auto inner_product = direction_.cast<Complex>().dot(E);
            std::abs(inner_product) > 1e-10 * E.norm()) {
            throw std::invalid_argument("E must be transverse to launch direction");
        }
    }

    void propagate(double distance, double refractiveIndex);

    // Call after propagating to a perfect electric conductor (PEC) surface.
    // Save the outgoing field at the same position and optical path as the hit.
    void reflect_pec(const Vec3 &unitNormal);
};

Vec3 reflected_direction(
    const Vec3 &incidentDirection,
    const Vec3 &unitNormal);

// Returns the entry distance into the reception sphere, or 0 if already inside.
// receiverPosition must be finite and radius must be finite and positive.
std::optional<double> intersect_receiver(
    const EMRay &ray,
    const Vec3 &receiverPosition,
    double radius);


enum class TraceStatus { Received, Escaped, ReflectionLimit, DistanceLimit, AmbiguousHit };

struct TraceOptions {
    int maxReflections_ = 8;
    double maxDistance_ = 20.0; // physical distance in metres, from this call's start
    double refractiveIndex_ = 1.0;
};

struct TraceResult {
    TraceStatus status_;
    int reflections_;
    double distance_;
};

// Trace one ray. After the last permitted reflection, still check the next
// segment for reception. Escaping rays are drawn out to maxDistance_.
// When launched inside/on the Rx sphere, reception starts after leaving it.
TraceResult trace_ray(EMRay &ray, const std::vector<Surface> &surfaces,
                      const Vec3 &receiverPosition, double receiverRadius,
                      const TraceOptions &options = {});


std::vector<Vec3> launch_directions(std::size_t count);

Vec3C launch_polarization(const Vec3 &unitDirection);
