#pragma once

#include "EMRay.hpp"

// Open an interactive 3D view of the saved paths and polarization ellipses.
void visualize_rays(const std::vector<EMRay> &rays,
                    const std::vector<Rectangle> &rectangles,
                    const Vec3 &receiverPosition, double receiverRadius);
