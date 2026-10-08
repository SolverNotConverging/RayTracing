#pragma once

#include "EMRay.hpp"

// Open an interactive 3D view of the saved paths and polarization ellipses.
// Surfaces must have been validated during scene setup.
void visualize_rays(const std::vector<EMRay> &rays,
                    const std::vector<Surface> &surfaces,
                    const Vec3 &receiverPosition, double receiverRadius);

// Geometry-only view for refined paths, before field reconstruction.
void visualize_paths(const std::vector<std::vector<Vec3>> &paths,
                     const std::vector<Surface> &surfaces,
                     const Vec3 &receiverPosition, double receiverRadius,
                     const std::vector<std::vector<Vec3>> &unresolvedPaths = {});
