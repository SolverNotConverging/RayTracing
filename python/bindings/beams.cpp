#include "bindings.hpp"
#include <RayTracing/BeamSolver.hpp>
#include <pybind11/functional.h>

void bind_beams(py::module_ &root) {
    namespace b = rt::beam;
    auto m=root.def_submodule("_beams", "Shared C++ Gaussian-beam numerics");
    py::class_<b::Material>(m,"Material")
        .def(py::init([](std::string kind,Complex epsilon,double mu,double conductivity,std::string name) {
            b::Material value{std::move(kind),epsilon,mu,conductivity,std::move(name)}; value.validate(); return value;
        }),py::arg("kind"),py::arg("epsilon_r")=Complex(1),py::arg("mu_r")=1.,py::arg("conductivity")=0.,py::arg("name")="")
        .def_readonly("kind",&b::Material::kind).def_readonly("epsilon_r",&b::Material::epsilon)
        .def_readonly("mu_r",&b::Material::mu).def_readonly("conductivity",&b::Material::conductivity)
        .def_readonly("name",&b::Material::name).def("permittivity",&b::Material::permittivity)
        .def("index",&b::Material::index).def("validate",&b::Material::validate);
    m.def("fresnel",&b::fresnel);
    py::class_<b::Box>(m,"Box").def(py::init<Vec3,Vec3,Eigen::Matrix3d>())
        .def_readwrite("center",&b::Box::center).def_readwrite("half_lengths",&b::Box::halfLengths)
        .def_readwrite("axes",&b::Box::axes).def("validate",&b::Box::validate);
    m.def("faces",&b::faces); m.def("contains",&b::contains); m.def("edge_distance",&b::edge_distance);
    m.def("curvature",&b::curvature); m.def("unit",&b::unit); m.def("frame",&b::frame); m.def("polarization",&b::polarization);
    py::class_<b::Object>(m,"Object").def_readonly("shape",&b::Object::shape)
        .def_readonly("material",&b::Object::material).def_readonly("volume",&b::Object::volume).def_readonly("terminate",&b::Object::terminate);
    py::class_<b::Event>(m,"Event").def(py::init<>())
        .def_readwrite("distance",&b::Event::distance).def_readwrite("position",&b::Event::position)
        .def_readwrite("normal",&b::Event::normal).def_readwrite("shape",&b::Event::shape)
        .def_readwrite("object_id",&b::Event::objectId).def_readwrite("before",&b::Event::before)
        .def_readwrite("after",&b::Event::after).def_readwrite("kind",&b::Event::kind);
    py::class_<b::Scene>(m,"Scene").def(py::init<b::Material,double>(),py::arg("background")=b::Material{},py::arg("tolerance")=1e-8)
        .def("add",&b::Scene::add,py::arg("shape"),py::arg("material"),py::arg("volume"),py::arg("terminate")="")
        .def("region_at",&b::Scene::region_at).def("material",&b::Scene::material)
        .def("starts_terminated",&b::Scene::starts_terminated).def("next_event",&b::Scene::next_event)
        .def("validate",&b::Scene::validate).def_property_readonly("background",&b::Scene::background)
        .def_property_readonly("tolerance",&b::Scene::tolerance).def_property_readonly("objects",&b::Scene::objects);
    py::class_<b::State>(m,"State").def(py::init<>())
        .def_readwrite("position",&b::State::position).def_readwrite("direction",&b::State::direction)
        .def_readwrite("basis",&b::State::basis).def_readwrite("curvature",&b::State::curvature)
        .def_readwrite("field",&b::State::field).def_readwrite("region",&b::State::region)
        .def_readwrite("power",&b::State::power).def_readwrite("distance",&b::State::distance)
        .def_readwrite("optical_path",&b::State::opticalPath).def_readwrite("history",&b::State::history)
        .def_readwrite("branch",&b::State::branch).def_readwrite("parent",&b::State::parent).def("validate",&b::State::validate);
    m.def("propagate",[](const b::State &s,double distance,const b::Material &mat,double frequency) {
        return b::propagate(s,distance,mat,frequency);
    });
    py::class_<b::Segment>(m,"Segment").def(py::init<b::State,double,b::Material,double>())
        .def_readwrite("state",&b::Segment::state).def_readwrite("length",&b::Segment::length)
        .def_readwrite("material",&b::Segment::material).def_readwrite("frequency",&b::Segment::frequency)
        .def("field_at",&b::Segment::field_at).def("radii",&b::Segment::radii,py::arg("distance")=0.);
    py::class_<b::GaussianBeam>(m,"GaussianBeam").def(py::init<>())
        .def_readwrite("waist_position",&b::GaussianBeam::waistPosition).def_readwrite("direction",&b::GaussianBeam::direction)
        .def_readwrite("waist_radius",&b::GaussianBeam::waistRadius).def_readwrite("power",&b::GaussianBeam::power)
        .def_readwrite("polarization",&b::GaussianBeam::polarization).def_readwrite("launch_distance",&b::GaussianBeam::launchDistance)
        .def_readwrite("phase",&b::GaussianBeam::phase).def_readwrite("transverse_axis",&b::GaussianBeam::transverseAxis)
        .def("launch",&b::GaussianBeam::launch).def("beamlets",&b::GaussianBeam::beamlets,
            py::arg("samples")=17,py::arg("radius_ratio")=.5,py::arg("extent")=4.);
    py::class_<b::PlaneWave>(m,"PlaneWave").def(py::init<>())
        .def_readwrite("position",&b::PlaneWave::position).def_readwrite("direction",&b::PlaneWave::direction)
        .def_readwrite("polarization",&b::PlaneWave::polarization).def_readwrite("amplitude",&b::PlaneWave::amplitude)
        .def("launch",&b::PlaneWave::launch);
    py::class_<b::ApertureSource>(m,"ApertureSource").def(py::init<>())
        .def_readwrite("position",&b::ApertureSource::position).def_readwrite("direction",&b::ApertureSource::direction)
        .def_readwrite("half_widths",&b::ApertureSource::halfWidths).def_readwrite("samples",&b::ApertureSource::samples)
        .def_readwrite("beamlet_radius",&b::ApertureSource::beamletRadius).def_readwrite("field",&b::ApertureSource::field)
        .def_readwrite("rcond",&b::ApertureSource::rcond).def_readwrite("frequency_hz",&b::ApertureSource::frequencyHz)
        .def_readwrite("vacuum_only",&b::ApertureSource::vacuumOnly).def("launch",&b::ApertureSource::launch)
        .def_static("from_antenna",&b::ApertureSource::from_antenna);
    py::class_<b::SolverConfig>(m,"SolverConfig").def(py::init<>())
        .def_readwrite("frequency_hz",&b::SolverConfig::frequencyHz).def_readwrite("max_interactions",&b::SolverConfig::maxInteractions)
        .def_readwrite("max_distance",&b::SolverConfig::maxDistance).def_readwrite("min_power_fraction",&b::SolverConfig::minPowerFraction)
        .def_readwrite("max_branches",&b::SolverConfig::maxBranches).def_readwrite("max_loss_ratio",&b::SolverConfig::maxLossRatio)
        .def("validate",&b::SolverConfig::validate);
    py::class_<b::Terminal>(m,"Terminal").def(py::init<>())
        .def_readwrite("branch",&b::Terminal::branch).def_readwrite("reason",&b::Terminal::reason)
        .def_readwrite("position",&b::Terminal::position).def_readwrite("power",&b::Terminal::power).def_readwrite("history",&b::Terminal::history);
    py::class_<b::Simulation>(m,"Simulation").def(py::init<>())
        .def_readwrite("config",&b::Simulation::config).def_readwrite("scene",&b::Simulation::scene)
        .def_readwrite("segments",&b::Simulation::segments).def_readwrite("terminals",&b::Simulation::terminals)
        .def_readwrite("warnings",&b::Simulation::warnings).def_readwrite("launched_power",&b::Simulation::launchedPower)
        .def_readwrite("absorbed_power",&b::Simulation::absorbedPower).def_readwrite("interface_flux_residual",&b::Simulation::interfaceFluxResidual)
        .def("field_at",&b::Simulation::field_at,py::call_guard<py::gil_scoped_release>())
        .def("fields_at",&b::Simulation::fields_at,py::call_guard<py::gil_scoped_release>())
        .def("power_balance",&b::Simulation::power_balance);
    m.def("interact",&b::interact);
    // Snapshot under the GIL, then release it for the numerical transport loop.
    // Python aperture callbacks are evaluated during launch, before this call.
    m.def("solve_states",[](b::Scene scene,std::vector<b::State> states,b::SolverConfig config) {
        py::gil_scoped_release release; return b::solve_states(scene,states,config);
    });
}
