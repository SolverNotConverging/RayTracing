#include <RayTracing/Solver.hpp>
#include <RayTracing/Constants.hpp>
#ifdef RAYTRACING_EXAMPLE_IO
#include <RayTracing/ResultIO.hpp>
#endif
#include "visualization.hpp"
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) try {
    rt::SolverConfig config;
    config.rayCount = 20000;
    std::filesystem::path output, plot, screenshot;
    bool interactive = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (++i == argc) throw std::invalid_argument("Missing value for " + arg);
            return argv[i];
        };
        if (arg == "--help") {
            std::cout << "point_ray_example [--rays N]\n"
                         "  [--output DIR]       HDF5/CSV (requires IO build)\n"
                         "  [--plot FILE.png]    Impulse magnitude/phase (requires Matplot++)\n"
                         "  [--view]            Interactive ray viewer (requires VTK)\n"
                         "  [--screenshot FILE.png] Offscreen image, or save before --view\n";
            return 0;
        } else if (arg == "--rays") {
            const auto count = value();
            if (count.empty() || count.find_first_not_of("0123456789") != std::string::npos)
                throw std::invalid_argument("Ray count must be a positive integer");
            config.rayCount = std::stoull(count);
            if (config.rayCount == 0) throw std::invalid_argument("Ray count must be positive");
        } else if (arg == "--output") output = value();
        else if (arg == "--plot") plot = value();
        else if (arg == "--screenshot") screenshot = value();
        else if (arg == "--view") interactive = true;
        else throw std::invalid_argument("Unknown option: " + arg);
    }
#ifndef RAYTRACING_EXAMPLE_IO
    if (!output.empty()) throw std::invalid_argument("Rebuild with RAYTRACING_BUILD_IO=ON for --output");
#endif
#ifndef RAYTRACING_EXAMPLE_MATPLOT
    if (!plot.empty()) throw std::invalid_argument("Rebuild with RAYTRACING_EXAMPLE_MATPLOT=ON for --plot");
#endif
#ifndef RAYTRACING_EXAMPLE_VTK
    if (interactive || !screenshot.empty())
        throw std::invalid_argument("Rebuild with RAYTRACING_EXAMPLE_VTK=ON for --view/--screenshot");
#endif
    const Vec3 x = Vec3::UnitX(), y = Vec3::UnitY(), z = Vec3::UnitZ();
    const auto pec = rt::Material::pec();
    rt::Scene scene;
    scene.add(Rectangle{{1.5, 0, 0}, y, z, 3, 1}, pec);
    scene.add(Rectangle{{-1.5, 0, 0}, y, z, 3, 1}, pec);
    scene.add(Rectangle{{0, -3, 0}, x, z, 3, 1}, pec);
    scene.add(Rectangle{{0, 3, 0}, x, z, 3, 1}, pec);
    scene.add(Sphere{{0.7, -1.2, 0}, 0.3}, pec);
    scene.add(Disk{{-0.7, -1.3, 0}, y, 0.35}, pec);
    scene.add(Cylinder{{0.6, 2, 0}, Vec3(0.2, 0, 1).normalized(), 0.25, 0.7, true}, pec);
    scene.add(Triangle{{-1, 0, 0.8}, {0, 0, 0.8}, {-0.5, 1, 0.8}}, pec);
    const rt::Antenna antenna = rt::ThinWireDipole{rt::constants::speedOfLight / config.frequencyHz / 2, 1};
    const rt::Transmitter tx{Vec3::Zero(), antenna};
    const rt::Receiver rx{{-0.5, 1.4, 0}, rt::Isotropic{}};
    const auto result = rt::solve(scene, tx, rx, config);
    std::cout << result.rays.size() << " paths, " << result.responseRayIndices.size()
              << " paths with fields, " << result.impulseResponse.taps_.size() << " impulse taps\n";
    std::cout << "delay_seconds,coefficient_real,coefficient_imag\n" << std::setprecision(17);
    for (const auto &tap : result.impulseResponse.taps_)
        std::cout << tap.delaySeconds_ << ',' << tap.coefficient_.real() << ',' << tap.coefficient_.imag() << '\n';
#ifdef RAYTRACING_EXAMPLE_IO
    if (!output.empty()) {
        std::filesystem::create_directories(output);
        rt::save_h5(result, output / "simulation.h5");
        rt::save_impulse_csv(result.impulseResponse, output / "impulse_response.csv");
    }
#endif
    for (const auto &file : {plot, screenshot})
        if (file.has_parent_path()) std::filesystem::create_directories(file.parent_path());
#ifdef RAYTRACING_EXAMPLE_MATPLOT
    if (!plot.empty()) plot_impulse(result.impulseResponse, plot);
#endif
#ifdef RAYTRACING_EXAMPLE_VTK
    if (interactive || !screenshot.empty()) view_rays(result, interactive, screenshot);
#endif
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
