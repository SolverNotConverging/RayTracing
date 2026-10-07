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
    constexpr double FREQUENCY = 77e9;
    constexpr double WAVELENGTH0 = C / FREQUENCY;
    const Vec3 initial_position = Vec3(0.0, 0.0, 0.0);
    const Vec3 initial_direction = Vec3{-1.0, 0.0, 0.0};
    const Vec3C E = Vec3{0.0, 1.0, 0.0};

    const Vec3 RxCenter(2.0, 0.0, 0.0);
    const double RxRadius{0.1};

    EMRay ray(initial_position, initial_direction, E, WAVELENGTH0);
    const auto planePoint = Vec3{4.0, 0.0, 0.0};
    const auto planeNormal = Vec3{-1.0, 0.0, 0.0};


    const auto receiverDistance = intersect_receiver(ray, RxCenter, RxRadius);
    const auto planeDistance = intersect_plane(ray, planePoint, planeNormal);

    // On a tie, the surface takes priority. Never detect Rx through a wall.
    if (receiverDistance && (!planeDistance || *receiverDistance < *planeDistance)) {
        ray.propagate(*receiverDistance, 1.0);
        std::cout << "Entered reception sphere at distance: " << *receiverDistance << " m\n";
        std::cout << "Entry position:\n" << ray.path_.back().position << '\n';
        std::cout << "Optical path to sphere entry: " << ray.path_.back().opticalPath << " m\n";
    } else if (planeDistance) {
        std::cout << "Hit PEC plane at distance: " << *planeDistance << " m\n";
        ray.propagate(*planeDistance, 1.0);
        ray.reflect_pec(planeNormal);
    } else {
        std::cout << "No Hit" << std::endl;
    }
}
