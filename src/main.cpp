#include <EMRay.hpp>
#include "Visualization.hpp"
#include <cmath>
#include <iostream>


Vec3 launch_direction(double theta, double phi) {
    double x = std::sin(theta) * std::cos(phi);
    double y = std::sin(theta) * std::sin(phi);
    double z = std::cos(theta);
    return Vec3(x, y, z);
}

Vec3C launch_E(double theta, double phi) {
    double Ex = std::cos(theta) * std::cos(phi);
    double Ey = std::cos(theta) * std::sin(phi);
    double Ez = -std::sin(theta);
    return Vec3(Ex, Ey, Ez);
}

int main() {
    std::vector<EMRay> rays;

    const Vec3 source_position = Vec3(0.0, 0.0, 0.0);

    constexpr int NUM_PHI = 5;
    constexpr int NUM_THETA = 5;
    constexpr int NUM_RAYS = NUM_PHI * NUM_THETA;

    constexpr double FREQUENCY = 77e9;
    constexpr double WAVELENGTH0 = C / FREQUENCY;


    for (int iTheta = 0; iTheta < NUM_THETA; ++iTheta) {
        const double theta =
                (static_cast<double>(iTheta) + 0.5) /
                static_cast<double>(NUM_THETA) *
                std::numbers::pi;

        for (int iPhi = 0; iPhi < NUM_PHI; ++iPhi) {
            const double phi =
                    static_cast<double>(iPhi) /
                    static_cast<double>(NUM_PHI) *
                    2.0 * std::numbers::pi;

            const Vec3 direction =
                    launch_direction(theta, phi);

            const Vec3C E =
                    launch_E(theta, phi);

            rays.emplace_back(
                source_position,
                direction,
                E,
                WAVELENGTH0
            );
        }
    }

    for (int i = 0; i < NUM_RAYS; ++i) {
        // Save several samples along each 2 m ray.
        for (int step = 0; step < 5; ++step) {
            propagate(rays[i], 0.4, 1.0);
        }
        std::cout << "ray index" << i + 1 << std::endl;
        std::cout << "position" << rays[i].path_.back().position << std::endl;
        std::cout << "direction" << rays[i].direction_ << std::endl;
    }
    visualize_rays(rays);
}
