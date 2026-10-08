#pragma once
#include <Eigen/Dense>
#include <complex>
#include <limits>
#include <optional>
#include <variant>
#include <vector>

using Vec3 = Eigen::Vector3d;
using Complex = std::complex<double>;
using Vec3C = Eigen::Vector3cd;

struct Rectangle {
    Vec3 center_;
    Vec3 u_; // Unit width direction
    Vec3 v_; // Unit height direction, perpendicular to u
    double halfWidth_;
    double halfHeight_;
};

struct Disk {
    Vec3 center_;
    Vec3 normal_; // Unit normal
    double radius_;
};

struct Sphere {
    Vec3 center_;
    double radius_;
};

struct Cylinder {
    Vec3 center_; // Midpoint of the axis
    Vec3 axis_; // Unit axis; ends are center_ +/- halfLength_ * axis_
    double radius_;
    double halfLength_;
    bool capped_ = true;
};

struct Triangle {
    Vec3 a_;
    Vec3 b_;
    Vec3 c_; // Vertex order determines the normal
};

using Surface = std::variant<Rectangle, Disk, Sphere, Cylinder, Triangle>;

struct GeometryHit {
    double distance_;
    Vec3 normal_; // Unit geometric normal, not flipped to face the ray
};

struct SurfaceHit {
    double distance_;
    Vec3 normal_;
    std::size_t surfaceIndex_;
    bool ambiguous_ = false; // Simultaneous hits with incompatible normals
    std::vector<std::size_t> coincidentSurfaceIndices_;
};

// Independent of the self-hit cutoff tMin. Distances are measured along a unit ray.
inline constexpr double simultaneousHitTolerance = 1e-8; // metres
inline constexpr double equivalentNormalTolerance = 1e-8; // unit-vector chord length

// Call once after scene construction, and again if geometry is changed.
// Throws for nonfinite geometry, invalid dimensions, or invalid axes/vertices.
void validate_surface(const Surface &surface);

// Shapes must already be validated; intersection routines do not recheck geometry.
// All surfaces are two-sided. Direction must be unit length. Return the nearest
// hit in (tMin, tMax]; inside a sphere/cylinder this can be the exit surface.
// Sphere/cylinder normals point outward. At a cylinder rim, caps win exact ties.
std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Rectangle &shape, double tMin = 1e-8,
                                     double tMax = std::numeric_limits<double>::infinity());

std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Disk &shape, double tMin = 1e-8,
                                     double tMax = std::numeric_limits<double>::infinity());

std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Sphere &shape, double tMin = 1e-8,
                                     double tMax = std::numeric_limits<double>::infinity());

std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Cylinder &shape, double tMin = 1e-8,
                                     double tMax = std::numeric_limits<double>::infinity());

std::optional<GeometryHit> intersect(const Vec3 &origin, const Vec3 &direction,
                                     const Triangle &shape, double tMin = 1e-8,
                                     double tMax = std::numeric_limits<double>::infinity());

std::optional<SurfaceHit> nearest_surface(const Vec3 &origin, const Vec3 &direction,
                                          const std::vector<Surface> &surfaces, double tMin = 1e-8,
                                          double tMax = std::numeric_limits<double>::infinity());
