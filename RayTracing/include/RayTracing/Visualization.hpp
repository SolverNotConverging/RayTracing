#pragma once

#include "EMRay.hpp"
#include "Solver.hpp"

// Open an interactive 3D view of the saved paths and polarization ellipses.
// Surfaces must have been validated during scene setup.
void visualize_rays(const std::vector<EMRay> &rays,
                    const std::vector<Surface> &surfaces,
                    const Vec3 &receiverPosition, double receiverRadius);

// Geometry-only view for refined paths, before field reconstruction.
void visualize_paths(
    const std::vector<std::vector<Vec3> > &paths,
    const std::vector<Surface> &surfaces, const Vec3 &receiverPosition,
    double receiverRadius,
    const std::vector<std::vector<Vec3> > &unresolvedPaths = {});

namespace rt {
    struct ViewOptions {
        bool interactive = true;
        std::filesystem::path screenshot;
        double patternScale =
                0.45; // Peak display radius (m), independent of physical fields.
        std::vector<std::size_t> rayIndices; // Empty selects every solved ray.
    };

    void visualize(const SimulationResult &result, const ViewOptions &options = {});

    void visualize(const SimulationResult &result, std::size_t rayIndex,
                   const ViewOptions &options = {});

    void visualize_h5(const std::filesystem::path &file,
                      const ViewOptions &options = {});

    void plot(const ImpulseResponse &response,
              const std::filesystem::path &png = "impulse_response.png",
              bool showWindow = true);

    void plot_csv(const std::filesystem::path &csv,
                  const std::filesystem::path &png = "impulse_response.png",
                  bool showWindow = true);
} // namespace rt
