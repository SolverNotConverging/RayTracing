#include <raytracing/geometry.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace rt {
    namespace {
        constexpr double angularTolerance = 1e-12;
        constexpr double edgeTolerance = 1e-9; // metres, planar face boundaries

        void require(bool condition, const char *message) {
            if (!condition) throw std::invalid_argument(message);
        }

        bool unit(const Vec3 &v) { return v.allFinite() && std::abs(v.norm() - 1.0) <= 1e-10; }
        bool positive(double x) { return std::isfinite(x) && x > 0.0; }

        void check(const Rectangle &s) {
            require(s.center.allFinite() && unit(s.u) && unit(s.v) && std::abs(s.u.dot(s.v)) <= 1e-10 &&
                    positive(s.halfWidth) && positive(s.halfHeight),
                    "Rectangle needs a finite centre, perpendicular unit axes and positive half sizes");
        }

        void check(const Disk &s) {
            require(s.center.allFinite() && unit(s.normal) && positive(s.radius),
                    "Disk needs a finite centre, unit normal and positive radius");
        }

        void check(const Sphere &s) {
            require(s.center.allFinite() && positive(s.radius), "Sphere needs a finite centre and positive radius");
        }

        void check(const Cylinder &s) {
            require(s.center.allFinite() && unit(s.axis) && positive(s.radius) && positive(s.halfLength),
                    "Cylinder needs a finite centre, unit axis and positive radius and half length");
        }

        void check(const Triangle &s) {
            const Vec3 ab = s.b - s.a, ac = s.c - s.a;
            require(s.a.allFinite() && s.b.allFinite() && s.c.allFinite() &&
                    ab.cross(ac).norm() > angularTolerance * ab.norm() * ac.norm(),
                    "Triangle vertices must be finite, distinct and non-collinear");
        }

        void check(const Box &s) {
            require(s.center.allFinite() && s.halfLengths.allFinite() && (s.halfLengths.array() > 0).all() &&
                    (s.axes.transpose() * s.axes - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() <= 1e-12 &&
                    s.axes.determinant() > 0,
                    "Box needs positive half lengths and right-handed orthonormal axes");
        }

        bool in_range(double t, double tMin, double tMax) { return t > tMin && t <= tMax; }

        std::optional<double> plane_distance(const Vec3 &origin, const Vec3 &direction, const Vec3 &point,
                                             const Vec3 &normal, double tMin, double tMax) {
            const double denominator = normal.dot(direction);
            if (std::abs(denominator) < angularTolerance) return std::nullopt;
            const double t = normal.dot(point - origin) / denominator;
            return in_range(t, tMin, tMax) ? std::optional<double>{t} : std::nullopt;
        }

        double sign(double x) { return (x > 0) - (x < 0); }
    } // namespace

    void validate_surface(const Surface &surface) {
        std::visit([](const auto &s) { check(s); }, surface);
    }

    void validate_shape(const Shape &shape) {
        std::visit([](const auto &s) { check(s); }, shape);
    }

    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Rectangle &s, double tMin,
                                 double tMax) {
        const Vec3 normal = s.u.cross(s.v).normalized();
        const auto t = plane_distance(origin, direction, s.center, normal, tMin, tMax);
        if (!t) return std::nullopt;
        const Vec3 offset = origin + *t * direction - s.center;
        if (std::abs(offset.dot(s.u)) > s.halfWidth + edgeTolerance ||
            std::abs(offset.dot(s.v)) > s.halfHeight + edgeTolerance)
            return std::nullopt;
        return Hit{*t, normal};
    }

    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Disk &s, double tMin,
                                 double tMax) {
        const auto t = plane_distance(origin, direction, s.center, s.normal, tMin, tMax);
        if (!t || (origin + *t * direction - s.center).norm() > s.radius + edgeTolerance) return std::nullopt;
        return Hit{*t, s.normal};
    }

    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Sphere &s, double tMin,
                                 double tMax) {
        const Vec3 offset = origin - s.center;
        const double middle = -offset.dot(direction);
        const Vec3 perpendicular = offset + middle * direction;
        const double halfChordSquared = s.radius * s.radius - perpendicular.squaredNorm();
        if (halfChordSquared < 0.0) return std::nullopt;
        const double halfChord = std::sqrt(halfChordSquared);
        for (double t: {middle - halfChord, middle + halfChord})
            if (in_range(t, tMin, tMax)) return Hit{t, (offset + t * direction).normalized()};
        return std::nullopt;
    }

    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Cylinder &s, double tMin,
                                 double tMax) {
        const Vec3 offset = origin - s.center;
        const double axialOrigin = offset.dot(s.axis);
        const double axialDirection = direction.dot(s.axis);
        const Vec3 radialOrigin = offset - axialOrigin * s.axis;
        const Vec3 radialDirection = direction - axialDirection * s.axis;
        std::optional<Hit> nearest;
        auto consider = [&](double t, const Vec3 &normal) {
            if (in_range(t, tMin, tMax) && (!nearest || t < nearest->distance)) nearest = Hit{t, normal};
        };

        // Caps first, so their normals win an exact rim tie.
        if (s.capped && std::abs(axialDirection) >= angularTolerance) {
            for (double side: {-1.0, 1.0}) {
                const double t = (side * s.halfLength - axialOrigin) / axialDirection;
                if ((radialOrigin + t * radialDirection).norm() <= s.radius + edgeTolerance) consider(t, side * s.axis);
            }
        }

        // Circle intersection in the plane perpendicular to the axis.
        const double a = radialDirection.squaredNorm();
        if (a > 0.0) {
            const double middle = -radialOrigin.dot(radialDirection) / a;
            const Vec3 closest = radialOrigin + middle * radialDirection;
            const double halfChordSquared = (s.radius * s.radius - closest.squaredNorm()) / a;
            if (halfChordSquared >= 0.0) {
                const double halfChord = std::sqrt(halfChordSquared);
                for (double t: {middle - halfChord, middle + halfChord})
                    if (std::abs(axialOrigin + t * axialDirection) <= s.halfLength)
                        consider(t, (radialOrigin + t * radialDirection).normalized());
            }
        }
        return nearest;
    }

    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Triangle &s, double tMin,
                                 double tMax) {
        const Vec3 ab = s.b - s.a, ac = s.c - s.a;
        const Vec3 areaNormal = ab.cross(ac);
        const Vec3 normal = areaNormal.normalized();
        const auto t = plane_distance(origin, direction, s.a, normal, tMin, tMax);
        if (!t) return std::nullopt;
        const Vec3 offset = origin + *t * direction - s.a;
        const double areaSquared = areaNormal.squaredNorm();
        const double beta = offset.cross(ac).dot(areaNormal) / areaSquared;
        const double gamma = ab.cross(offset).dot(areaNormal) / areaSquared;
        constexpr double barycentricTolerance = 1e-10;
        if (beta < -barycentricTolerance || gamma < -barycentricTolerance || beta + gamma > 1.0 + barycentricTolerance)
            return std::nullopt;
        return Hit{*t, normal};
    }

    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Surface &shape, double tMin,
                                 double tMax) {
        return std::visit([&](const auto &s) { return intersect(origin, direction, s, tMin, tMax); }, shape);
    }

    std::optional<SurfaceHit> nearest_surface(const Vec3 &origin, const Vec3 &direction,
                                              const std::vector<Surface> &surfaces, double tMin, double tMax) {
        std::vector<SurfaceHit> hits;
        for (std::size_t i = 0; i < surfaces.size(); ++i)
            if (const auto hit = intersect(origin, direction, surfaces[i], tMin, tMax))
                hits.push_back({hit->distance, hit->normal, i});
        if (hits.empty()) return std::nullopt;
        SurfaceHit nearest = *std::min_element(hits.begin(), hits.end(), [](const auto &a, const auto &b) {
            return a.distance < b.distance;
        });
        // Group against the actual minimum, not transitively.
        std::vector<Vec3> normals;
        for (const auto &hit: hits) {
            if (hit.distance - nearest.distance > simultaneousHitTolerance) continue;
            nearest.coincidentSurfaces.push_back(hit.surfaceIndex);
            for (const Vec3 &normal: normals)
                // Opposite normals describe the same two-sided plane.
                if (std::min((normal - hit.normal).norm(), (normal + hit.normal).norm()) > equivalentNormalTolerance)
                    nearest.ambiguous = true;
            normals.push_back(hit.normal);
        }
        return nearest;
    }

    std::vector<Surface> faces(const Shape &shape) {
        return std::visit([](const auto &s) -> std::vector<Surface> {
            if constexpr (std::is_same_v<std::decay_t<decltype(s)>, Box>) {
                std::vector<Surface> result;
                for (int i = 0; i < 3; ++i) {
                    const int j = (i + 1) % 3, k = (i + 2) % 3;
                    for (double side: {-1.0, 1.0})
                        result.push_back(Rectangle{s.center + side * s.halfLengths[i] * s.axes.col(i), s.axes.col(j),
                                                   side * s.axes.col(k), s.halfLengths[j], s.halfLengths[k]});
                }
                return result;
            } else return {s};
        }, shape);
    }

    bool contains(const Shape &shape, const Vec3 &p) {
        return std::visit([&](const auto &s) -> bool {
            using T = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<T, Box>)
                return ((s.axes.transpose() * (p - s.center)).cwiseAbs().array() < s.halfLengths.array()).all();
            else if constexpr (std::is_same_v<T, Sphere>) return (p - s.center).norm() < s.radius;
            else if constexpr (std::is_same_v<T, Cylinder>) {
                const Vec3 v = p - s.center;
                const double a = v.dot(s.axis);
                return std::abs(a) < s.halfLength && (v - a * s.axis).norm() < s.radius;
            } else return false;
        }, shape);
    }

    double edge_distance(const Surface &shape, const Vec3 &p) {
        return std::visit([&](const auto &s) -> double {
            using T = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<T, Rectangle>) {
                const Vec3 v = p - s.center;
                return std::min(s.halfWidth - std::abs(v.dot(s.u)), s.halfHeight - std::abs(v.dot(s.v)));
            } else if constexpr (std::is_same_v<T, Disk>) return s.radius - (p - s.center).norm();
            else if constexpr (std::is_same_v<T, Cylinder>) {
                const Vec3 v = p - s.center;
                const double axial = v.dot(s.axis);
                return std::hypot(s.halfLength - std::abs(axial), s.radius - (v - axial * s.axis).norm());
            } else if constexpr (std::is_same_v<T, Triangle>) {
                const auto edge = [&](const Vec3 &a, const Vec3 &b) { return (b - a).cross(p - a).norm() / (b - a).norm(); };
                return std::min({edge(s.a, s.b), edge(s.b, s.c), edge(s.c, s.a)});
            } else return infinity;
        }, shape);
    }

    Eigen::Matrix3d curvature(const Surface &shape, const Vec3 &p, const Vec3 &normal) {
        return std::visit([&](const auto &s) -> Eigen::Matrix3d {
            using T = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<T, Sphere>) {
                const Vec3 outward = (p - s.center) / s.radius;
                return sign(outward.dot(normal)) * (Eigen::Matrix3d::Identity() - outward * outward.transpose()) /
                       s.radius;
            } else if constexpr (std::is_same_v<T, Cylinder>) {
                const Vec3 v = p - s.center, radial = v - v.dot(s.axis) * s.axis;
                if (radial.norm() > 0 && std::abs(normal.dot(s.axis)) < 0.5) {
                    const Vec3 outward = radial.normalized();
                    return sign(outward.dot(normal)) *
                           (Eigen::Matrix3d::Identity() - s.axis * s.axis.transpose() - outward * outward.transpose()) /
                           s.radius;
                }
            }
            return Eigen::Matrix3d::Zero();
        }, shape);
    }
} // namespace rt
