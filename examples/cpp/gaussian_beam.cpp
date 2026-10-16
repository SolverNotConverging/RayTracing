#include <RayTracing/BeamSolver.hpp>
#include <iomanip>
#include <iostream>

// Numerical example only: applications can plot these complex samples externally.
int main() try {
    namespace b=rt::beam;
    b::Scene scene;
    scene.add(b::Box{{0,0,.5065},{2,2,.0065}}, b::Material{"dielectric",{4,.02},1,0,"glass"},true);
    b::GaussianBeam source;
    source.waistRadius.setConstant(.1);
    b::SolverConfig config;
    config.frequencyHz=10e9; config.maxDistance=4; config.maxInteractions=32; config.minPowerFraction=0;
    const auto result=b::solve(scene,b::Source{source},config);
    std::cout << "x,y,z,Ex_real,Ex_imag,Ey_real,Ey_imag,Ez_real,Ez_imag\n" << std::setprecision(17);
    for (double z : {.2,.7,1.,1.5}) {
        const Vec3 p(.03,.01,z); const auto e=result.field_at(p);
        std::cout << p.x() << ',' << p.y() << ',' << p.z();
        for (int i=0;i<3;++i) std::cout << ',' << e[i].real() << ',' << e[i].imag();
        std::cout << '\n';
    }
    return 0;
} catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
