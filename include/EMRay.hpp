#pragma once
#include <Eigen/Dense>
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
        direction_(direction.normalized()),
        wavelength0_(wavelength0) {

        if (E.norm() < 1e-10) {
            throw std::invalid_argument{"E cannot be zero"};
        }
        if (const auto inner_product = direction.cast<Complex>().dot(E); std::abs(inner_product) > 1e-10) {
            throw std::invalid_argument("E must be transverse to launch direction");
        }
    }
};

std::optional<double> intersect_plane(
    const EMRay &ray,
    const Vec3 &planePoint,
    const Vec3 &planeNormal);

void propagate(EMRay &ray,
               double distance,
               double refractiveIndex);
