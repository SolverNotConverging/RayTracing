// Multipath impulse response between two antennas in a corner reflector.
#include <raytracing/solver.hpp>
#include <iostream>

int main() {
    rt::Scene scene;
    scene.add(rt::Rectangle{{2, 0, 0}, {0, 1, 0}, {0, 0, 1}, 3, 1}, rt::Material::pec());
    scene.add(rt::Rectangle{{0, 2, 0}, {1, 0, 0}, {0, 0, 1}, 3, 1}, rt::Material::conductor(5.8e7));

    rt::Transmitter tx{{0, 0, 0}, rt::ThinWireDipole{rt::speedOfLight / 77e9 / 2}};
    rt::Receiver rx{{1, 0.5, 0}, rt::Isotropic{}};
    const auto result = rt::solve(scene, tx, rx);

    for (const auto &tap: result.impulseResponse.taps)
        std::cout << tap.delay * 1e9 << " ns  " << tap.coefficient << '\n';
}
