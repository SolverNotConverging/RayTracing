#include <RayTracing/BeamSolver.hpp>
#include <RayTracing/Constants.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace b=rt::beam;
namespace {
constexpr double pi=rt::constants::pi, c0=rt::constants::speedOfLight, eta=rt::constants::freeSpaceImpedance;
void check(bool ok,const char *message) { if (!ok) throw std::runtime_error(message); }
void near(Complex actual,Complex expected,double tolerance,const char *message) {
    check(std::abs(actual-expected)<=tolerance*std::max(1.,std::abs(expected)),message);
}
template<class F> void rejects(F action) { bool rejected=false; try { action(); } catch (const std::invalid_argument &) { rejected=true; } check(rejected,"Invalid native input accepted"); }
b::Material dielectric(Complex epsilon) { return {"dielectric",epsilon,1,0,"dielectric"}; }
b::Material pec() { return {"pec",1,1,0,"PEC"}; }
}
int main() try {
    const double frequency=10e9, w=.1, k=2*pi*frequency/c0, zr=k*w*w/2, e0=std::sqrt(4*eta/(pi*w*w));
    b::SolverConfig config; config.frequencyHz=frequency; config.maxDistance=4;
    b::GaussianBeam source; source.waistPosition={0,0,1}; source.waistRadius.setConstant(w); source.launchDistance=-1;
    const auto free=b::solve(b::Scene{},b::Source{source},config);
    for (const double z : {.1,.6,1.,1.7,3.}) {
        const double s=z-1;
        const auto expected=e0/(1.+Complex(0,s/zr))*std::exp(Complex(0,k*s));
        near(free.field_at({0,0,z})[0],expected,2e-12,"Vacuum Gaussian amplitude/Gouy phase changed");
    }
    const auto state=free.segments[0].state;
    auto many=state; for (int i=0;i<100;++i) many=b::propagate(many,.02,{},frequency);
    check((many.field-b::propagate(state,2,{},frequency).field).norm()<1e-9,"Propagation lost semigroup phase");
    check(free.field_at({0,0,4}).norm()==0,"Segment endpoint must be half open");

    // Independently sum the scalar Fabry-Perot series for a lossy slab.
    for (const Complex epsilon : {Complex(4),Complex(4,.02)}) {
        b::Scene slab; slab.add(b::Box{{0,0,.5065},{2,2,.0065}},dielectric(epsilon),true);
        b::PlaneWave wave; auto cfg=config; cfg.maxInteractions=32; cfg.minPowerFraction=0;
        const auto result=b::solve(slab,b::Source{wave},cfg);
        const Complex n=std::sqrt(epsilon), p=std::exp(Complex(0,k)*n*.013), r=(n-1.)/(n+1.);
        const Complex expected=2./(1.+n)*2.*n/(n+1.)*p/(1.-r*r*p*p)*std::exp(Complex(0,k*(1-.013)));
        near(result.field_at({0,0,1})[0],expected,2e-12,"Coherent slab phase/refraction incorrect");
        const auto budget=result.power_balance(); double accounted=0;
        for (const auto &[name,power] : budget) if (name!="launched") accounted+=power;
        near(accounted,result.launchedPower,1e-12,"Branch power budget does not close");
    }

    // Concave PEC sphere focuses a finite beam; the focus must retain Gouy phase.
    b::Scene mirror; mirror.add(Sphere{Vec3::Zero(),2},pec(),false);
    b::GaussianBeam beam; beam.waistRadius.setConstant(w);
    auto cfg=config; cfg.maxDistance=3.99; cfg.maxInteractions=1;
    const auto focused=b::solve(mirror,b::Source{beam},cfg);
    const auto &reflected=focused.segments.at(1);
    const Complex q0(0,-zr), qi=q0+2., qo=1./(1./qi-1.);
    for (const double distance : {0.,.5,1.,1.9}) {
        const auto actual=b::propagate(reflected.state,distance,{},frequency);
        const auto expected=-e0*q0/qi*qo/(qo+distance)*std::exp(Complex(0,k*(2+distance)));
        near(actual.field[0],expected,2e-11,"Curved-mirror field/Gouy phase incorrect");
        const auto radii=reflected.radii(distance);
        near(radii[0],std::sqrt(2/(k*(1./(qo+distance)).imag())),2e-12,"Focused beam width incorrect");
    }

    b::Scene regions;
    const auto outer=regions.add(b::Box{Vec3::Zero(),Vec3::Constant(2)},dielectric(4),true);
    const auto inner=regions.add(Sphere{Vec3::Zero(),1},dielectric(9),true);
    check(regions.region_at(Vec3::Zero())==static_cast<int>(inner) && regions.region_at({1.5,0,0})==static_cast<int>(outer),"Region insertion priority changed");
    regions.add(b::Box{{0,0,0},{3,3,3}},std::nullopt,true,"entry");
    check(b::solve(regions,b::Source{b::PlaneWave{}},config).terminals[0].reason=="terminated","Volume termination ignored");

    const auto direct=b::solve(b::Scene{},b::Source{beam},config);
    std::vector<b::Source> beamlets;
    for (const auto &s : beam.beamlets(25)) beamlets.emplace_back(s);
    cfg=config; cfg.minPowerFraction=0;
    const auto split=b::solve(b::Scene{},beamlets,cfg);
    const Vec3 probe(.03,-.02,.7);
    check((split.field_at(probe)-direct.field_at(probe)).norm()/direct.field_at(probe).norm()<2e-6,"Native beamlet reconstruction failed");
    b::ApertureSource aperture; aperture.field=[w](const Vec3 &p) -> Vec3C { return Vec3C(std::exp(-(p.x()*p.x()+p.y()*p.y())/(w*w)),0,0); };
    const auto fitted=b::solve(b::Scene{},b::Source{aperture},cfg);
    near(fitted.field_at(Vec3::Zero())[0],1,1e-8,"Native aperture least-squares fit failed");
    const auto adapted=b::ApertureSource::from_antenna(rt::Antenna{rt::Isotropic{}},{0,0,-10},Vec3::Zero(),Vec3::UnitZ(),{.2,.2},5,.1,frequency);
    check(b::solve(b::Scene{},b::Source{adapted},config).field_at({0,0,.1}).allFinite(),"Native antenna aperture adapter failed");

    auto bad=config; bad.maxLossRatio=.1; rejects([&] { b::solve(b::Scene{},b::Source{beam},bad); });
    auto invalid=beam; invalid.waistRadius[0]=0; rejects([&] { invalid.launch(b::Scene{},frequency); });
    rejects([&] { b::Scene{}.add(Rectangle{Vec3::Zero(),Vec3::UnitX(),Vec3::UnitY(),1,1},dielectric(4),true); });
    std::cout << "Native Gaussian propagation, coherent slabs, focus phase, regions, beamlets and apertures passed\n";
    return 0;
} catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
