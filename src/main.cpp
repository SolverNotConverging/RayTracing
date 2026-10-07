#include <EMRay.hpp>
#include "Visualization.hpp"
#include <cmath>
#include <iostream>


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

int main() {
    std::vector<EMRay> rays;
    constexpr double FREQUENCY = 77e9;
    constexpr double WAVELENGTH0 = C / FREQUENCY;
    const Vec3 initial_position = Vec3(0.0, 0.0, 0.5);
    const Vec3 initial_direction = Vec3{3.0 / 5.0, 0.0, 4.0 / 5.0};
    const Vec3C E = Vec3{0.0, 1.0, 0.0};

    const Vec3 planePoint(0.0, 0.0, 2.0);
    const Vec3 planeNormal(0.0, 0.0, 1.0);

    EMRay ray(initial_position, initial_direction, E, WAVELENGTH0);

    auto distance = intersect_plane(ray, planePoint, planeNormal);
    if (distance.has_value()) {
        std::cout << "Distance: " << *distance << std::endl;
        propagate(ray, *distance, 1.0);
    } else {
        std::cout << "No forward intersection" << std::endl;
    }
}
