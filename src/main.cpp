#include <EMRay.hpp>
#include "Visualization.hpp"
#include "TraceConfig.hpp"
#include <iostream>
#include <utility>

int main(int argc, char *argv[]) try {
    if (argc > 2) {
        std::cerr << "Usage: RayTracing [trace_options.json]\n";
        return 1;
    }
    const std::filesystem::path configPath =
            argc == 2 ? argv[1] : RAYTRACING_DEFAULT_CONFIG;
    const TraceOptions options = load_trace_options(configPath);
    constexpr double FREQUENCY = 77e9;
    constexpr double WAVELENGTH0 = C / FREQUENCY;

    const Vec3 RxCenter{-0.5, 1.4, 0.0};
    constexpr double RxRadius = 0.12;
    const std::vector<Surface> surfaces{
        Rectangle{Vec3{1.5, 0.0, 0.0}, Vec3::UnitY(), Vec3::UnitZ(), 3.0, 1.0},
        Rectangle{Vec3{-1.5, 0.0, 0.0}, Vec3::UnitY(), Vec3::UnitZ(), 3.0, 1.0},
        Rectangle{Vec3{0, -3.0, 0.0}, Vec3::UnitX(), Vec3::UnitZ(), 3.0, 1.0},
        Rectangle{Vec3{0, 3.0, 0.0}, Vec3::UnitX(), Vec3::UnitZ(), 3.0, 1.0},
        // Additional analytical PEC shapes; curved surfaces are meshed only in VTK.
        Sphere{Vec3{0.7, -1.2, 0.0}, 0.3},
        Disk{Vec3{-0.7, -1.3, 0.0}, Vec3::UnitY(), 0.35},
        Cylinder{Vec3{0.6, 2.0, 0.0}, Vec3{0.2, 0.0, 1.0}.normalized(), 0.25, 0.7, true},
        Triangle{Vec3{-1.0, 0.0, 0.8}, Vec3{0.0, 0.0, 0.8}, Vec3{-0.5, 1.0, 0.8}},
    };

    // Geometry is constant during tracing: validate each shape once at setup.
    for (const Surface &surface: surfaces) {
        validate_surface(surface);
    }

    constexpr std::size_t NUM_RAYS = 2000;

    std::vector<EMRay> receivedRays;

    for (const Vec3 &direction: launch_directions(NUM_RAYS)) {
        EMRay ray(
            Vec3::Zero(),
            direction,
            launch_polarization(direction),
            WAVELENGTH0);

        const TraceResult result =
                trace_ray(ray, surfaces, RxCenter, RxRadius, options);

        if (result.status_ == TraceStatus::Received) {
            receivedRays.push_back(std::move(ray));
        }
    }

    std::cout << "Received candidates: "
            << receivedRays.size() << '\n' << std::flush;

    visualize_rays(receivedRays, surfaces, RxCenter, RxRadius);
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
