#include "Tracing.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
    void check(bool condition, const char *message) {
        if (!condition) throw std::runtime_error(message);
    }

    void hit_is(const std::optional<GeometryHit> &hit, double distance, const Vec3 &normal) {
        check(hit.has_value(), "Expected an intersection");
        check(std::abs(hit->distance_ - distance) < 1e-8, "Incorrect hit distance");
        check((hit->normal_ - normal).norm() < 1e-9, "Incorrect local normal");
        check(std::abs(hit->normal_.norm() - 1.0) < 1e-10, "Normal is not unit length");
    }

    template<class Action>
    void rejects(Action action) {
        bool rejected = false;
        try { action(); } catch (const std::invalid_argument &) { rejected = true; }
        check(rejected, "Expected invalid geometry or query to be rejected");
    }
}

int main() try {
    const Vec3 X = Vec3::UnitX(), Y = Vec3::UnitY(), Z = Vec3::UnitZ();
    const Rectangle rectangle{Vec3(0, 0, 2), X, Y, 1, 2};
    hit_is(intersect(Vec3::Zero(), Z, rectangle), 2, Z);
    hit_is(intersect(Vec3(1, 2, 0), Z, rectangle), 2, Z);
    check(!intersect(Vec3(1.1, 0, 0), Z, rectangle), "Rectangle bounds ignored");
    check(!intersect(Vec3::Zero(), X, rectangle), "Parallel rectangle ray accepted");
    hit_is(intersect(Vec3(0, 0, 3), -Z, rectangle), 1, Z);
    const Disk disk{Vec3(0, 0, 2), Z, 1};
    hit_is(intersect(Vec3(1, 0, 0), Z, disk), 2, Z);
    check(!intersect(Vec3(1, 1, 0), Z, disk), "Disk used square instead of circular bounds");
    check(!intersect(Vec3(0, 0, 2), X, disk), "Coplanar disk ray accepted");

    const Sphere sphere{Vec3(2, 0, 0), 1};
    hit_is(intersect(Vec3::Zero(), X, sphere), 1, -X);
    hit_is(intersect(sphere.center_, X, sphere), 1, X);
    hit_is(intersect(Vec3(1, 0, 0), X, sphere), 2, X);
    check(!intersect(Vec3(1, 0, 0), -X, sphere), "Outward self-hit accepted");
    hit_is(intersect(Vec3(0, 1, 0), X, sphere), 2, Y);
    check(!intersect(Vec3(0, 1.1, 0), X, sphere), "Sphere miss accepted");
    check(!intersect(Vec3::Zero(), -X, sphere), "Behind-ray sphere accepted");
    hit_is(intersect(Vec3::Zero(), X, Sphere{Vec3(1e8, 0, 0), 0.1}), 1e8 - 0.1, -X);
    check(!intersect(Vec3::Zero(), X, sphere, 0, 0.9), "tMax ignored");
    hit_is(intersect(Vec3::Zero(), X, sphere, 0, 1), 1, -X);
    hit_is(intersect(Vec3::Zero(), X, sphere, 1, 3), 3, X);
    check(!intersect(Vec3::Zero(), X, sphere, 1, 1), "tMin must be exclusive");

    const Cylinder cylinder{Vec3::Zero(), Z, 1, 1, true};
    hit_is(intersect(Vec3(2, 0, 0), -X, cylinder), 1, X);
    hit_is(intersect(Vec3::Zero(), X, cylinder), 1, X);
    hit_is(intersect(Vec3(0, 0, 3), -Z, cylinder), 2, Z);
    hit_is(intersect(Vec3(0, 0, -3), Z, cylinder), 2, -Z);
    hit_is(intersect(Vec3::Zero(), Z, cylinder), 1, Z);
    check(!intersect(Vec3(2, 0, 3), -Z, cylinder), "Cap radius ignored");
    hit_is(intersect(Vec3(2, 1, 0), -X, cylinder), 2, Y);
    const Cylinder open{Vec3::Zero(), Z, 1, 1, false};
    check(!intersect(Vec3::Zero(), Z, open), "Open cylinder unexpectedly has caps");
    check(!intersect(Vec3(2, 0, 2), -X, open), "Cylinder height ignored");
    hit_is(intersect(Vec3(2, 0, 3), Vec3(-1, 0, -1).normalized(), open), 3 * std::sqrt(2.0), -X);
    hit_is(intersect(Vec3(2, 0, 3), Vec3(-1, 0, -1).normalized(), cylinder), 2 * std::sqrt(2.0), Z);
    hit_is(intersect(Vec3(2, 0, 2), Vec3(-1, 0, -1).normalized(), cylinder), std::sqrt(2.0), Z);

    const Triangle triangle{Vec3(0, 0, 2), Vec3(2, 0, 2), Vec3(0, 2, 2)};
    hit_is(intersect(Vec3(0.5, 0.5, 0), Z, triangle), 2, Z);
    hit_is(intersect(Vec3(1, 1, 0), Z, triangle), 2, Z);
    hit_is(intersect(Vec3::Zero(), Z, triangle), 2, Z);
    check(!intersect(Vec3(1.5, 1.5, 0), Z, triangle), "Triangle bounds ignored");
    hit_is(intersect(Vec3(0.5, 0.5, 3), -Z, triangle), 1, Z);
    hit_is(intersect(Vec3(0.5, 0.5, 0), Z, Triangle{triangle.a_, triangle.c_, triangle.b_}), 2, -Z);

    // Rotation and translation must preserve intersection distance and rotate normals.
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(0.7, Vec3(1, 2, 3).normalized()).toRotationMatrix();
    const Vec3 shift(4, -2, 7);
    hit_is(intersect(rotation * Vec3(2, 0, 0) + shift, rotation * -X,
                     Cylinder{shift, rotation * Z, 1, 1, true}), 1, rotation * X);
    hit_is(intersect(shift, rotation * Z,
                     Disk{rotation * disk.center_ + shift, rotation * Z, 1}), 2, rotation * Z);
    hit_is(intersect(shift, rotation * Z,
                     Rectangle{rotation * rectangle.center_ + shift, rotation * X, rotation * Y, 1, 2}), 2,
           rotation * Z);
    hit_is(intersect(shift, rotation * X, Sphere{rotation * sphere.center_ + shift, 1}), 1, rotation * -X);
    hit_is(intersect(rotation * Vec3(0.5, 0.5, 0) + shift, rotation * Z,
                     Triangle{
                         rotation * triangle.a_ + shift, rotation * triangle.b_ + shift, rotation * triangle.c_ + shift
                     }), 2, rotation * Z);

    const std::vector<Surface> mixed{
        Rectangle{Vec3(5, 0, 0), Y, Z, 2, 2},
        Sphere{Vec3(3, 0, 0), 0.5},
        Disk{Vec3(1, 0, 0), X, 1},
        Cylinder{Vec3(7, 0, 0), Z, 1, 1, true},
        Triangle{Vec3(9, -1, -1), Vec3(9, 1, -1), Vec3(9, 0, 1)}
    };
    for (const Surface &surface: mixed) validate_surface(surface);
    const auto nearest = nearest_surface(Vec3::Zero(), X, mixed);
    check(nearest && nearest->surfaceIndex_ == 2 && nearest->distance_ == 1, "Mixed nearest selection failed");
    check(!nearest_surface(Vec3::Zero(), X, mixed, 0, 0.9), "Nearest query ignored range");
    check(!nearest_surface(Vec3::Zero(), X, {}), "Empty scene returned a hit");

    rejects([&] { validate_surface(Rectangle{Vec3::Zero(), 2 * X, 0.5 * Y, 1, 1}); });
    rejects([&] { validate_surface(Rectangle{Vec3::Zero(), X, X, 1, 1}); });
    rejects([&] { validate_surface(Disk{Vec3::Zero(), Z, 0}); });
    rejects([&] { validate_surface(Sphere{Vec3::Zero(), -1}); });
    rejects([&] { validate_surface(Cylinder{Vec3::Zero(), Vec3::Zero(), 1, 1}); });
    rejects([&] { validate_surface(Triangle{Vec3::Zero(), X, 2 * X}); });
    rejects([&] { intersect(Vec3::Zero(), 2 * X, sphere); });
    rejects([&] { intersect(Vec3::Zero(), X, sphere, 2, 1); });

    // A sphere's local normal must give the same reflection law as a planar surface.
    TracedRay ray(Vec3::Zero(), X);
    const auto result = trace_ray(ray, {Sphere{Vec3(2, 0, 0), 0.5}}, Vec3(-1, 0, 0), 0.1, {1, 10});
    check(result.status_ == TraceStatus::Received && result.reflections_ == 1, "Curved reflection failed to reach Rx");
    check(std::abs(result.distance_ - 3.9) < 1e-10 && std::abs(ray.distance_ - 3.9) < 1e-10,
          "Reflected distance incorrect");
    check(ray.vertices_.size() == 3 && (ray.vertices_[1] - Vec3(1.5, 0, 0)).norm() < 1e-12,
          "Hit vertex history incorrect");
    check(ray.reflections_.size() == 1 && ray.reflections_[0].vertexIndex_ == 1, "Reflection event incorrect");
    const auto oblique = intersect(Vec3(0, 0.5, 0), X, sphere);
    TracedRay angled(Vec3(0, 0.5, 0), X);
    angled.advance(oblique->distance_);
    angled.reflect(oblique->normal_, 0);
    check(std::abs(angled.direction_.dot(oblique->normal_) + X.dot(oblique->normal_)) < 1e-10, "Reflection law failed");

    // An exact sharp corner must stop without choosing an arbitrary normal.
    const Surface cornerX = Rectangle{Vec3(1, 0, 0), Y, Z, 1, 1};
    const Surface cornerY = Rectangle{Vec3(0, 1, 0), X, Z, 1, 1};
    validate_surface(cornerX);
    validate_surface(cornerY);
    const Vec3 cornerDirection = Vec3(1, 1, 0).normalized();
    for (bool reversed : {false, true}) {
        const std::vector<Surface> corner = reversed ? std::vector<Surface>{cornerY, cornerX}
                                                   : std::vector<Surface>{cornerX, cornerY};
        TracedRay cornerRay(Vec3::Zero(), cornerDirection);
        // Keep Rx off the ray paths, and leave enough distance/reflection budget.
        const auto cornerResult = trace_ray(cornerRay, corner, Vec3(0, 0, 10), 0.1, {8, 20});
        check(cornerResult.status_ == TraceStatus::AmbiguousHit, "Exact corner not flagged");
        check(cornerResult.reflections_ == 0 && cornerRay.reflections_.empty(),
              "Ambiguous corner produced a reflection");
        check((cornerRay.vertices_[1] - Vec3(1, 1, 0)).norm() < 1e-12, "Ray missed corner");
        check(cornerRay.vertices_.size() == 2 && (cornerRay.direction_ - cornerDirection).norm() < 1e-12,
              "Ambiguous ray must stop with incident direction and sample");
        check(std::abs(cornerResult.distance_ - std::sqrt(2.0)) < 1e-12, "Corner stop distance incorrect");
        std::cout << "Exact corner (" << (reversed ? "y wall first" : "x wall first")
                  << "): AmbiguousHit, reflections=" << cornerResult.reflections_
                  << ", direction=" << cornerRay.direction_.transpose() << '\n';
    }
    // Nearby rays encounter two distinct points and reflect from both faces.
    for (double offset : {-1e-4, 1e-4}) {
        const Vec3 direction = Vec3(1, 1 + offset, 0).normalized();
        TracedRay nearby(Vec3::Zero(), direction);
        const auto nearbyResult = trace_ray(nearby, {cornerX, cornerY}, Vec3(0, 0, 10), 0.1, {8, 20});
        check(nearbyResult.status_ == TraceStatus::Escaped && nearbyResult.reflections_ == 2,
              "Near-corner ray should reflect twice");
        check((nearby.direction_ + direction).norm() < 1e-12, "Two perpendicular reflections should reverse this ray");
        std::cout << "Near corner (offset=" << offset << "): Escaped, reflections=2\n";
    }

    // A near tie uses its own tolerance, independent of the query's tMin.
    const auto nearTie = nearest_surface(Vec3::Zero(), Vec3(1, 1 + 1e-9, 0).normalized(),
                                         {cornerX, cornerY}, 1e-12);
    check(nearTie && nearTie->ambiguous_, "Near-simultaneous corner not detected");
    TracedRay beforeCorner(Vec3::Zero(), cornerDirection);
    check(trace_ray(beforeCorner, {cornerX, cornerY}, Vec3(0.5, 0.5, 0), 0.01).status_ == TraceStatus::Received,
          "Corner beyond receiver prevented reception");
    TracedRay shortCorner(Vec3::Zero(), cornerDirection);
    check(trace_ray(shortCorner, {cornerX, cornerY}, Vec3(0, 0, 10), 0.1, {8, 1}).status_ == TraceStatus::DistanceLimit,
          "Corner beyond distance limit terminated trace early");

    // Two coplanar triangles with opposite winding meet at the ray hit.
    const std::vector<Surface> seam{
        Triangle{Vec3(1, -1, -1), Vec3(1, 1, -1), Vec3(1, 1, 1)},
        Triangle{Vec3(1, -1, -1), Vec3(1, -1, 1), Vec3(1, 1, 1)}};
    for (const auto &surface : seam) validate_surface(surface);
    const auto seamHit = nearest_surface(Vec3::Zero(), X, seam);
    check(seamHit && !seamHit->ambiguous_ && seamHit->coincidentSurfaceIndices_.size() == 2,
          "Coplanar seam incorrectly ambiguous");
    TracedRay seamRay(Vec3::Zero(), X);
    const auto seamResult = trace_ray(seamRay, seam, -X, 0.1);
    check(seamResult.status_ == TraceStatus::Received && seamResult.reflections_ == 1,
          "Coplanar seam failed to reflect normally");

    std::cout << "All surface and mixed-scene tracing checks passed\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
