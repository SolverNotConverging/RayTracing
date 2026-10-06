#pragma once

#include "crt/ray.hpp"

#include <cmath>
#include <optional>

namespace crt {

    struct Triangle {
        Vec3 v0;
        Vec3 v1;
        Vec3 v2;
    };

    struct Hit {
        double distance;
        Vec3 position;
        Vec3 normal;

        double u;
        double v;
    };

    inline std::optional<Hit>
    intersect(const Ray& ray, const Triangle& triangle)
    {
        constexpr double eps = 1e-10;

        const Vec3 edge1 = triangle.v1 - triangle.v0;
        const Vec3 edge2 = triangle.v2 - triangle.v0;

        const Vec3 pvec = ray.direction.cross(edge2);

        const double det = edge1.dot(pvec);

        if (std::abs(det) < eps) {
            return std::nullopt;
        }

        const double invDet = 1.0 / det;

        const Vec3 tvec = ray.origin - triangle.v0;

        const double u = tvec.dot(pvec) * invDet;

        if (u < 0.0 || u > 1.0) {
            return std::nullopt;
        }

        const Vec3 qvec = tvec.cross(edge1);

        const double v = ray.direction.dot(qvec) * invDet;

        if (v < 0.0 || u + v > 1.0) {
            return std::nullopt;
        }

        const double t = edge2.dot(qvec) * invDet;

        if (t <= eps) {
            return std::nullopt;
        }

        Vec3 normal = edge1.cross(edge2).normalized();

        return Hit{
            .distance = t,
            .position = ray.at(t),
            .normal = normal,
            .u = u,
            .v = v
        };
    }

} // namespace crt