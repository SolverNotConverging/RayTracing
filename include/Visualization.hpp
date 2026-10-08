#pragma once

#include "EMRay.hpp"

// Open an interactive 3D view of the saved paths and polarization ellipses.
// Surfaces must have been validated during scene setup.
void visualize_rays(const std::vector<EMRay> &rays,
                    const std::vector<Surface> &surfaces,
                    const Vec3 &receiverPosition, double receiverRadius);
