#include <raytracing/beam.hpp>
#include <raytracing/core.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace b=rt::beam;
using namespace rt;

namespace {
constexpr double c0=speedOfLight, eta=freeSpaceImpedance;
void check(bool ok,const char *message) { if (!ok) throw std::runtime_error(message); }
void near(Complex actual,Complex expected,double tolerance,const char *message) {
    check(std::abs(actual-expected)<=tolerance*std::max(1.,std::abs(expected)),message);
}
template<class F> void rejects(F action) { bool rejected=false; try { action(); } catch (const std::invalid_argument &) { rejected=true; } check(rejected,"Invalid native input accepted"); }
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
        const auto response=free.impulse_response({0,0,z});
        check(response.taps.size()==1,"Vacuum beam must have one impulse tap");
        near(response.taps[0].delay*c0,z,1e-12,"Impulse delay must start at the launch plane");
        near(frequency_response(response),expected,2e-12,"Impulse sum must reconstruct the carrier field");
    }
    const auto state=free.segments[0].state;
    auto many=state; for (int i=0;i<100;++i) many=b::propagate(many,.02,{},frequency);
    check((many.field-b::propagate(state,2,{},frequency).field).norm()<1e-9,"Propagation lost semigroup phase");
    check(free.field_at({0,0,4}).norm()==0,"Segment endpoint must be half open");

    // Independently sum the scalar Fabry-Perot series for a lossy slab.
    for (const Complex epsilon : {Complex(4),Complex(4,.02)}) {
        b::Scene slab; slab.add_volume(Box{{0,0,.5065},{2,2,.0065}},Material::dielectric(epsilon));
        b::PlaneWave wave; auto cfg=config; cfg.maxInteractions=32; cfg.minPowerFraction=0;
        const auto result=b::solve(slab,b::Source{wave},cfg);
        const Complex n=std::sqrt(epsilon), p=std::exp(Complex(0,k)*n*.013), r=(n-1.)/(n+1.);
        const Complex expected=2./(1.+n)*2.*n/(n+1.)*p/(1.-r*r*p*p)*std::exp(Complex(0,k*(1-.013)));
        near(result.field_at({0,0,1})[0],expected,2e-12,"Coherent slab phase/refraction incorrect");
        const auto response=result.impulse_response({0,0,1});
        check(response.taps.size()>1,"Slab impulse response must resolve internal echoes");
        near(response.taps[0].delay*c0,1-.013+n.real()*.013,1e-12,"Slab optical delay incorrect");
        near(frequency_response(response),expected,2e-12,"Slab taps must retain coherent phase");
        const auto budget=result.power_balance(); double accounted=0;
        for (const auto &[name,power] : budget) if (name!="launched") accounted+=power;
        near(accounted,result.launchedPower,1e-12,"Branch power budget does not close");
    }

    // Concave PEC sphere focuses a finite beam; the focus must retain Gouy phase.
    b::Scene mirror; mirror.add_surface(Sphere{Vec3::Zero(),2},Material::pec());
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
    const auto outer=regions.add_volume(Box{Vec3::Zero(),Vec3::Constant(2)},Material::dielectric(4));
    const auto inner=regions.add_volume(Sphere{Vec3::Zero(),1},Material::dielectric(9));
    check(regions.region_at(Vec3::Zero())==static_cast<int>(inner) && regions.region_at({1.5,0,0})==static_cast<int>(outer),"Region insertion priority changed");
    regions.add_termination(Box{{0,0,0},{3,3,3}},true,b::Termination::Entry);
    check(b::solve(regions,b::Source{b::PlaneWave{}},config).terminals[0].reason==b::TerminalReason::Terminated,"Volume termination ignored");

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

    auto invalid=beam; invalid.waistRadius[0]=0; rejects([&] { invalid.launch(b::Scene{},frequency); });
    rejects([&] { b::Scene{}.add_volume(Rectangle{Vec3::Zero(),Vec3::UnitX(),Vec3::UnitY(),1,1},Material::dielectric(4)); });
    std::cout << "Native Gaussian propagation, coherent slabs, focus phase, regions, beamlets and apertures passed\n";
    return 0;
} catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
