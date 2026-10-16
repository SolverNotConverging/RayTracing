// Gaussian beam through a glass slab; sample the field behind it.
#include <raytracing/beam.hpp>
#include <iostream>

int main() {
    rt::beam::Scene scene;
    scene.add_volume(rt::Box{{0, 0, 0.5}, {1, 1, 0.01}}, rt::Material::dielectric({4, 0.02}));

    rt::beam::GaussianBeam beam;
    beam.waistRadius.setConstant(0.1);
    rt::beam::SolverConfig config;
    config.frequencyHz = 10e9;
    const auto result = rt::beam::solve(scene, beam, config);

    std::cout << "E(0, 0, 1) = " << result.field_at({0, 0, 1}).transpose() << '\n';
}
