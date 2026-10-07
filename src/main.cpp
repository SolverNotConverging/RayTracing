#include <EMRay.hpp>
#include "Visualization.hpp"
#include "TraceConfig.hpp"
#include <iostream>
#include <utility>


Vec3 launch_direction(double theta, double phi) {
    double x = std::sin(theta) * std::cos(phi);
    double y = std::sin(theta) * std::sin(phi);
    double z = std::cos(theta);
    return Vec3{x, y, z};
}

Vec3C launch_E(double theta, double phi) {
    double Ex = std::cos(theta) * std::cos(phi);
    double Ey = std::cos(theta) * std::sin(phi);
    double Ez = -std::sin(theta);
    return Vec3{Ex, Ey, Ez};
}

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
    const std::vector<Rectangle> rectangles{
        {Vec3{1.5, 0.0, 0.0}, Vec3::UnitY(), Vec3::UnitZ(), 3.0, 1.0},
        {Vec3{-1.5, 0.0, 0.0}, Vec3::UnitY(), Vec3::UnitZ(), 3.0, 1.0},
        {Vec3{0, -3.0, 0.0}, Vec3::UnitX(), Vec3::UnitZ(), 3.0, 1.0},
        {Vec3{0, 3.0, 0.0}, Vec3::UnitX(), Vec3::UnitZ(), 3.0, 1.0},
    };

    constexpr std::size_t NUM_RAYS = 2000;

    std::vector<EMRay> receivedRays;

    for (const Vec3 &direction: launch_directions(NUM_RAYS)) {
        EMRay ray(
            Vec3::Zero(),
            direction,
            launch_polarization(direction),
            WAVELENGTH0);

        const TraceResult result =
                trace_ray(ray, rectangles, RxCenter, RxRadius, options);

        if (result.status == TraceStatus::Received) {
            receivedRays.push_back(std::move(ray));
        }
    }

    std::cout << "Received candidates: "
            << receivedRays.size() << '\n';

    visualize_rays(receivedRays, rectangles, RxCenter, RxRadius);
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
