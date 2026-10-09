#include <RayTracing/ResultIO.hpp>
#include <RayTracing/Solver.hpp>
#ifdef RAYTRACING_HAS_VISUALIZATION
#include <RayTracing/Visualization.hpp>
#endif
#include <iostream>

int main(int argc, char **argv) try {
    bool interactive = true, isotropic = false;
    std::filesystem::path output = "results", input, csv, screenshot;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "--no-gui")
            interactive = false;
        else if (flag == "--isotropic")
            isotropic = true;
        else if (flag == "--output" && i + 1 < argc)
            output = argv[++i];
        else if (flag == "--load" && i + 1 < argc)
            input = argv[++i];
        else if (flag == "--plot-csv" && i + 1 < argc)
            csv = argv[++i];
        else if (flag == "--screenshot" && i + 1 < argc)
            screenshot = argv[++i];
        else
            throw std::invalid_argument(
                "Usage: RayTracing [--no-gui] [--isotropic] [--output dir] [--load "
                "simulation.h5] [--plot-csv response.csv] [--screenshot scene.png]");
    }
    if (!csv.empty()) {
#ifdef RAYTRACING_HAS_VISUALIZATION
        rt::plot_csv(csv, output / "impulse_response.png", interactive);
        return 0;
#else
        throw std::runtime_error("Visualization was disabled at build time");
#endif
    }

    rt::SolverConfig config;
    config.frequencyHz = 77e9;
    config.rayCount = 2000;
    config.maxReflections = 8;
    config.maxDistance = 20;
    rt::Scene scene;
    scene.add(
        rt::Rectangle{{1.5, 0, 0}, rt::Vec3::UnitY(), rt::Vec3::UnitZ(), 3, 1});
    scene.add(
        rt::Rectangle{{-1.5, 0, 0}, rt::Vec3::UnitY(), rt::Vec3::UnitZ(), 3, 1});
    scene.add(
        rt::Rectangle{{0, -3, 0}, rt::Vec3::UnitX(), rt::Vec3::UnitZ(), 3, 1});
    scene.add(
        rt::Rectangle{{0, 3, 0}, rt::Vec3::UnitX(), rt::Vec3::UnitZ(), 3, 1});
    scene.add(rt::Sphere{{0.7, -1.2, 0}, 0.3});
    scene.add(rt::Disk{{-0.7, -1.3, 0}, rt::Vec3::UnitY(), 0.35});
    scene.add(rt::Cylinder{
        {0.6, 2, 0}, rt::Vec3(0.2, 0, 1).normalized(), 0.25, 0.7, true
    });
    scene.add(rt::Triangle{{-1, 0, 0.8}, {0, 0, 0.8}, {-0.5, 1, 0.8}});
    rt::Transmitter tx{{0, 0, 0}, rt::Isotropic{}};
    rt::Receiver rx{{-0.5, 1.4, 0}, rt::Isotropic{}};
    if (input.empty() && !isotropic) {
        rt::FarfieldImportOptions import;
        import.inputConvention = rt::PhasorConvention::PositiveTime;
        tx.antenna = rt::load_farfield(RAYTRACING_EXAMPLE_PATTERN, import);
        std::cout << "CST source: " << tx.antenna.pattern->frequenciesHz.size()
                << " frequencies; " << tx.antenna.pattern->phiDegrees.size()
                << " x " << tx.antenna.pattern->thetaDegrees.size()
                << " samples\n";
    }
    const auto result =
            input.empty() ? rt::solve(scene, tx, rx, config) : rt::load_h5(input);
    rt::save_h5(result, output / "simulation.h5");
    rt::save_impulse_csv(result.impulseResponse, output / "impulse_response.csv");
    std::cout << "Launched: " << result.launchedRays
            << ", refined paths: " << result.rays.size()
            << ", impulse taps: " << result.impulseResponse.taps_.size()
            << '\n';
#ifdef RAYTRACING_HAS_VISUALIZATION
    rt::plot(result.impulseResponse, output / "impulse_response.png",
             interactive);
    if (interactive || !screenshot.empty())
        rt::visualize(result,
                      {.interactive = interactive, .screenshot = screenshot});
#endif
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
