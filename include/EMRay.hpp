#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>
#include <vector>
#include <optional>

using Vec3 = Eigen::Vector3d;
using Complex = std::complex<double>;
using Vec3C = Eigen::Vector3cd;

constexpr double PI = std::numbers::pi;
constexpr double C = 3e8;

struct RaySample {
    Vec3 position;
    Vec3C E;
    double opticalPath;
};

struct EMRay {
    std::vector<RaySample> path_;
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

// planeNormal must be a unit vector; returns a forward distance in metres.
std::optional<double> intersect_plane(
    const EMRay &ray,
    const Vec3 &planePoint,
    const Vec3 &planeNormal);

Vec3 reflected_direction(
    const Vec3 &incidentDirection,
    const Vec3 &unitNormal);
