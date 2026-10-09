#include "Surface.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace {
    constexpr double angularTolerance = 1e-12;
    constexpr double edgeTolerance = 1e-9; // metres, for planar face boundaries

    void require_point(const Vec3 &point) {
        if (!point.allFinite()) throw std::invalid_argument("Surface coordinates must be finite");
    }

    void require_unit(const Vec3 &axis) {
        if (!axis.allFinite() || std::abs(axis.norm() - 1.0) > 1e-10)
            throw std::invalid_argument("Surface axes and ray directions must be unit vectors");
    }

    void require_positive(double value) {
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument("Surface dimensions must be finite and positive");
    }

    void validate(const Rectangle &s) {
        require_point(s.center_);
        require_unit(s.u_);
        require_unit(s.v_);
        require_positive(s.halfWidth_);
        require_positive(s.halfHeight_);
        if (std::abs(s.u_.dot(s.v_)) > 1e-10)
            throw std::invalid_argument("Rectangle axes must be perpendicular");
    }

    void validate(const Disk &s) {
        require_point(s.center_);
        require_unit(s.normal_);
        require_positive(s.radius_);
    }

    void validate(const Sphere &s) {
        require_point(s.center_);
        require_positive(s.radius_);
    }

    void validate(const Cylinder &s) {
        require_point(s.center_);
        require_unit(s.axis_);
        require_positive(s.radius_);
        require_positive(s.halfLength_);
    }

    void validate(const Triangle &s) {
        require_point(s.a_);
        require_point(s.b_);
        require_point(s.c_);
        const Vec3 ab = s.b_ - s.a_, ac = s.c_ - s.a_;
        if (ab.norm() == 0.0 || ac.norm() == 0.0 ||
            ab.cross(ac).norm() <= angularTolerance * ab.norm() * ac.norm())
            throw std::invalid_argument("Triangle vertices must be distinct and non-collinear");
    }

    void validate_query(const Vec3 &origin, const Vec3 &direction, double tMin, double tMax) {
        require_point(origin);
        require_unit(direction);
        if (!std::isfinite(tMin) || tMin < 0.0 || std::isnan(tMax) || tMax < tMin)
            throw std::invalid_argument("Intersection interval must satisfy 0 <= tMin <= tMax");
    }

    bool in_range(double t, double tMin, double tMax) {
        return std::isfinite(t) && t > tMin && t <= tMax;
    }

    std::optional<double> plane_distance(const Vec3 &origin, const Vec3 &direction,
                                         const Vec3 &point, const Vec3 &normal, double tMin, double tMax) {
        const double denominator = normal.dot(direction);
        if (std::abs(denominator) < angularTolerance) return std::nullopt;
        const double t = normal.dot(point - origin) / denominator;
        return in_range(t, tMin, tMax) ? std::optional<double>{t} : std::nullopt;
    }
}

void validate_surface(const Surface &surface) {
    std::visit([](const auto &shape) { validate(shape); }, surface);
}

std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Rectangle &s, double tMin, double tMax) {
    validate_query(origin, direction, tMin, tMax);
    const Vec3 normal = s.u_.cross(s.v_).normalized();
    const auto t = plane_distance(origin, direction, s.center_, normal, tMin, tMax);
    if (!t) return std::nullopt;
    const Vec3 offset = origin + *t * direction - s.center_;
    if (std::abs(offset.dot(s.u_)) > s.halfWidth_ + edgeTolerance ||
        std::abs(offset.dot(s.v_)) > s.halfHeight_ + edgeTolerance)
        return std::nullopt;
    return GeometryHit{*t, normal};
}

std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Disk &s, double tMin, double tMax) {
    validate_query(origin, direction, tMin, tMax);
    const auto t = plane_distance(origin, direction, s.center_, s.normal_, tMin, tMax);
    if (!t || (origin + *t * direction - s.center_).norm() > s.radius_ + edgeTolerance)
        return std::nullopt;
    return GeometryHit{*t, s.normal_};
}

std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Sphere &s, double tMin, double tMax) {
    validate_query(origin, direction, tMin, tMax);
    const Vec3 offset = origin - s.center_;
    const double middle = -offset.dot(direction);
    const Vec3 perpendicular = offset + middle * direction;
    const double halfChordSquared = s.radius_ * s.radius_ - perpendicular.squaredNorm();
    if (halfChordSquared < 0.0) return std::nullopt;
    const double halfChord = std::sqrt(halfChordSquared);
    for (double t: {middle - halfChord, middle + halfChord}) {
        if (in_range(t, tMin, tMax))
            return GeometryHit{t, (offset + t * direction).normalized()};
    }
    return std::nullopt;
}

std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Cylinder &s, double tMin, double tMax) {
    validate_query(origin, direction, tMin, tMax);
    const Vec3 offset = origin - s.center_;
    const double axialOrigin = offset.dot(s.axis_);
    const double axialDirection = direction.dot(s.axis_);
    const Vec3 radialOrigin = offset - axialOrigin * s.axis_;
    const Vec3 radialDirection = direction - axialDirection * s.axis_;
    std::optional<GeometryHit> nearest;
    auto consider = [&](double t, const Vec3 &normal) {
        if (in_range(t, tMin, tMax) && (!nearest || t < nearest->distance_))
            nearest = GeometryHit{t, normal};
    };

    // Test caps first so that their normals win at an exact rim tie.
    if (s.capped_ && std::abs(axialDirection) >= angularTolerance) {
        for (double sign: {-1.0, 1.0}) {
            const double t = (sign * s.halfLength_ - axialOrigin) / axialDirection;
            if (in_range(t, tMin, tMax) &&
                (radialOrigin + t * radialDirection).norm() <= s.radius_ + edgeTolerance)
                consider(t, sign * s.axis_);
        }
    }

    // Solve the circle intersection in the plane perpendicular to the axis.
    const double a = radialDirection.squaredNorm();
    if (a > 0.0) {
        const double middle = -radialOrigin.dot(radialDirection) / a;
        const Vec3 closest = radialOrigin + middle * radialDirection;
        const double halfChordSquared = (s.radius_ * s.radius_ - closest.squaredNorm()) / a;
        if (halfChordSquared >= 0.0) {
            const double halfChord = std::sqrt(halfChordSquared);
            for (double t: {middle - halfChord, middle + halfChord}) {
                if (in_range(t, tMin, tMax) &&
                    std::abs(axialOrigin + t * axialDirection) <= s.halfLength_)
                    consider(t, (radialOrigin + t * radialDirection).normalized());
            }
        }
    }
    return nearest;
}

std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Triangle &s, double tMin, double tMax) {
    validate_query(origin, direction, tMin, tMax);
    const Vec3 ab = s.b_ - s.a_, ac = s.c_ - s.a_;
    const Vec3 areaNormal = ab.cross(ac);
    const Vec3 normal = areaNormal.normalized();
    const auto t = plane_distance(origin, direction, s.a_, normal, tMin, tMax);
    if (!t) return std::nullopt;
    const Vec3 offset = origin + *t * direction - s.a_;
    // Barycentric coordinates from oriented triangle areas.
    const double areaSquared = areaNormal.squaredNorm();
    const double beta = offset.cross(ac).dot(areaNormal) / areaSquared;
    const double gamma = ab.cross(offset).dot(areaNormal) / areaSquared;
    constexpr double barycentricTolerance = 1e-10;
    if (beta < -barycentricTolerance || gamma < -barycentricTolerance ||
        beta + gamma > 1.0 + barycentricTolerance)
        return std::nullopt;
    return GeometryHit{*t, normal};
}

std::optional<SurfaceHit> nearest_surface(const Vec3 &origin, const Vec3 &direction,
                                          const std::vector<Surface> &surfaces, double tMin, double tMax) {
    validate_query(origin, direction, tMin, tMax);
    std::vector<SurfaceHit> hits;
    for (std::size_t i = 0; i < surfaces.size(); ++i) {
        const auto hit = std::visit([&](const auto &shape) {
            return intersect(origin, direction, shape, tMin, tMax);
        }, surfaces[i]);
        if (hit) hits.push_back({hit->distance_, hit->normal_, i});
    }
    if (hits.empty()) return std::nullopt;
    SurfaceHit nearest = *std::min_element(hits.begin(), hits.end(),
                                           [](const SurfaceHit &a, const SurfaceHit &b) {
                                               return a.distance_ < b.distance_;
                                           });
    // Group against the actual minimum, not transitively against the previous hit.
    std::vector<Vec3> normals;
    for (const auto &hit: hits) {
        if (hit.distance_ - nearest.distance_ > simultaneousHitTolerance) continue;
        nearest.coincidentSurfaceIndices_.push_back(hit.surfaceIndex_);
        for (const Vec3 &normal: normals) {
            // Opposite normals describe the same two-sided reflecting plane.
            if (std::min((normal - hit.normal_).norm(), (normal + hit.normal_).norm()) > equivalentNormalTolerance)
                nearest.ambiguous_ = true;
        }
        normals.push_back(hit.normal_);
    }
    return nearest;
}
