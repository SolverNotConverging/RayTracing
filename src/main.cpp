#include <EMRay.hpp>
#include "Visualization.hpp"
#include "TraceConfig.hpp"
#include "PathRefinement.hpp"
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

    const Vec3 transmitterPosition = Vec3::Zero();
    RefinementOptions refinementOptions;
    refinementOptions.maxPathDistance_ = options.maxDistance_;
    std::vector<RefinementResult> refinedPaths;
    std::vector<EMRay> receivedCandidates; // Keep original sphere-entry geometry and reflection history.

    for (const Vec3 &direction: launch_directions(NUM_RAYS)) {
        EMRay ray(
            transmitterPosition,
            direction,
            launch_polarization(direction),
            WAVELENGTH0);

        const TraceResult result =
                trace_ray(ray, surfaces, RxCenter, RxRadius, options);

        if (result.status_ == TraceStatus::Received) {
            std::vector<std::size_t> sequence;
            for (const auto &event : ray.reflections_) sequence.push_back(event.surfaceIndex_);
            refinedPaths.push_back(refine_path(transmitterPosition, direction, RxCenter,
                                              surfaces, sequence, refinementOptions));
            receivedCandidates.push_back(std::move(ray));
        }
    }

    std::size_t converged = 0;
    std::cout << "Received candidates: " << refinedPaths.size() << '\n';
    for (std::size_t i = 0; i < refinedPaths.size(); ++i) {
        const auto &path = refinedPaths[i];
        if (path.status_ == RefinementStatus::Converged) ++converged;
        std::cout << "Candidate " << i << ": " << refinement_status_name(path.status_);
        if (path.geometry_.receiver_)
            std::cout << ", miss = " << path.geometry_.receiver_->missDistance_ << " m";
        std::cout << '\n';
    }
    const auto uniqueIndices = deduplicate_paths(refinedPaths);
    std::cout << "Converged: " << converged << ", unique paths: " << uniqueIndices.size()
              << '\n' << std::flush;

    // Display geometry only. Fields from the coarse rays do not belong to these paths.
    std::vector<std::vector<Vec3>> paths;
    std::vector<std::vector<Vec3>> unresolvedPaths;
    for (std::size_t i = 0; i < refinedPaths.size(); ++i) {
        if (refinedPaths[i].status_ != RefinementStatus::UnresolvedCorner) continue;
        std::vector<Vec3> points;
        for (const auto &sample : receivedCandidates[i].path_) points.push_back(sample.position_);
        unresolvedPaths.push_back(std::move(points));
    }
    std::cout << "Unresolved corner candidates retained: " << unresolvedPaths.size() << '\n' << std::flush;
    for (std::size_t index : uniqueIndices) {
        const auto &path = refinedPaths[index];
        std::vector<Vec3> points{path.transmitterPosition_};
        for (const auto &hit : path.geometry_.reflections_) points.push_back(hit.position_);
        points.push_back(path.geometry_.receiver_->closestPoint_);
        paths.push_back(std::move(points));
    }
    visualize_paths(paths, surfaces, RxCenter, RxRadius, unresolvedPaths);
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
