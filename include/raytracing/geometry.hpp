#pragma once
#include "core.hpp"
#include <limits>
#include <optional>
#include <variant>
#include <vector>

namespace rt {
    struct Rectangle {
        Vec3 center, u, v; // u and v: perpendicular unit axes; normal is u x v
        double halfWidth, halfHeight;
    };

    struct Disk {
        Vec3 center, normal; // unit normal
        double radius;
    };

    struct Sphere {
        Vec3 center;
        double radius;
    };

    struct Cylinder {
        Vec3 center, axis; // axis midpoint and unit axis; ends at center +/- halfLength * axis
        double radius, halfLength;
        bool capped = true;
    };

    struct Triangle {
        Vec3 a, b, c; // vertex order determines the normal
    };

    // Closed volume made of six rectangles; axes columns are a right-handed frame.
    struct Box {
        Vec3 center = Vec3::Zero(), halfLengths = Vec3::Ones();
        Eigen::Matrix3d axes = Eigen::Matrix3d::Identity();
    };

    using Surface = std::variant<Rectangle, Disk, Sphere, Cylinder, Triangle>;
    using Shape = std::variant<Box, Rectangle, Disk, Sphere, Cylinder, Triangle>;

    struct Hit {
        double distance;
        Vec3 normal; // unit geometric normal, not flipped to face the ray
    };

    struct SurfaceHit {
        double distance;
        Vec3 normal;
        std::size_t surfaceIndex;
        bool ambiguous = false; // simultaneous hits with incompatible normals
        std::vector<std::size_t> coincidentSurfaces;
    };

    inline constexpr double simultaneousHitTolerance = 1e-8; // metres
    inline constexpr double equivalentNormalTolerance = 1e-8; // unit-vector chord length
    inline constexpr double infinity = std::numeric_limits<double>::infinity();

    // Throw std::invalid_argument for nonfinite or degenerate geometry. Scenes call
    // these when shapes are added; the intersection routines assume valid shapes.
    void validate_surface(const Surface &surface);
    void validate_shape(const Shape &shape);

    // Nearest hit in (tMin, tMax] along a unit direction. Surfaces are two-sided;
    // sphere/cylinder normals point outward. At a cylinder rim, caps win ties.
    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Rectangle &shape,
                                 double tMin = 1e-8, double tMax = infinity);
    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Disk &shape,
                                 double tMin = 1e-8, double tMax = infinity);
    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Sphere &shape,
                                 double tMin = 1e-8, double tMax = infinity);
    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Cylinder &shape,
                                 double tMin = 1e-8, double tMax = infinity);
    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Triangle &shape,
                                 double tMin = 1e-8, double tMax = infinity);
    std::optional<Hit> intersect(const Vec3 &origin, const Vec3 &direction, const Surface &shape,
                                 double tMin = 1e-8, double tMax = infinity);

    std::optional<SurfaceHit> nearest_surface(const Vec3 &origin, const Vec3 &direction,
                                              const std::vector<Surface> &surfaces,
                                              double tMin = 1e-8, double tMax = infinity);

    // Boxes become six rectangles; other shapes are returned unchanged.
    std::vector<Surface> faces(const Shape &shape);
    // Strict interior test for boxes, spheres and capped cylinders; false for sheets.
    bool contains(const Shape &shape, const Vec3 &point);
    // Distance from a point on the surface to its nearest edge or rim (infinite for spheres).
    double edge_distance(const Surface &shape, const Vec3 &point);
    // Shape operator dn/dx at a surface point, for the given (either-sided) unit normal.
    Eigen::Matrix3d curvature(const Surface &shape, const Vec3 &point, const Vec3 &normal);
} // namespace rt
