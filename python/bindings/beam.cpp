// Gaussian-beam solver, exposed as raytracing.beam.
#include "bindings.hpp"

using namespace rt::beam;
using rt::Antenna, rt::Complex, rt::Material, rt::Vec3, rt::Vec3C;

namespace {
    using Points = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
    using Fields = Eigen::Matrix<Complex, Eigen::Dynamic, 3, Eigen::RowMajor>;

    // Accept a scalar or a pair for per-axis values.
    Eigen::Vector2d pair(const py::object &value) {
        if (py::isinstance<py::float_>(value) || py::isinstance<py::int_>(value))
            return Eigen::Vector2d::Constant(value.cast<double>());
        return value.cast<Eigen::Vector2d>();
    }

    std::vector<Source> sources(const py::object &value) {
        if (py::isinstance<py::list>(value) || py::isinstance<py::tuple>(value))
            return value.cast<std::vector<Source>>();
        return {value.cast<Source>()};
    }
} // namespace

void bind_beam(py::module_ &root) {
    auto m = root.def_submodule("beam", "Paraxial Gaussian-beam transport through material regions");

    m.def("frame", &frame, "direction"_a);
    py::enum_<Termination>(m, "Termination")
            .value("NONE", Termination::None).value("ENTRY", Termination::Entry).value("EXIT", Termination::Exit);
    py::enum_<EventKind>(m, "EventKind")
            .value("INTERFACE", EventKind::Interface).value("TERMINATE", EventKind::Terminate)
            .value("AMBIGUOUS", EventKind::Ambiguous);
    py::enum_<InteractionKind>(m, "InteractionKind")
            .value("REFLECT", InteractionKind::Reflect).value("TRANSMIT", InteractionKind::Transmit);
    py::enum_<TerminalReason>(m, "TerminalReason")
            .value("TERMINATED", TerminalReason::Terminated).value("POWER_CUTOFF", TerminalReason::PowerCutoff)
            .value("DISTANCE_LIMIT", TerminalReason::DistanceLimit)
            .value("INTERACTION_LIMIT", TerminalReason::InteractionLimit)
            .value("BRANCH_LIMIT", TerminalReason::BranchLimit).value("AMBIGUOUS", TerminalReason::Ambiguous)
            .value("GRAZING", TerminalReason::Grazing).value("LOSSY_CRITICAL", TerminalReason::LossyCritical);

    py::class_<Object>(m, "Object")
            .def_readonly("shape", &Object::shape).def_readonly("material", &Object::material)
            .def_readonly("volume", &Object::volume).def_readonly("termination", &Object::termination);
    py::class_<Event>(m, "Event")
            .def_readonly("distance", &Event::distance).def_readonly("position", &Event::position)
            .def_readonly("normal", &Event::normal).def_readonly("object", &Event::object)
            .def_readonly("before", &Event::before).def_readonly("after", &Event::after)
            .def_readonly("kind", &Event::kind);
    py::class_<Scene>(m, "Scene", "Material regions; the last added volume containing a point wins.")
            .def(py::init<Material, double>(), "background"_a = Material{}, "tolerance"_a = 1e-8)
            .def("add_volume", &Scene::add_volume, "shape"_a, "material"_a)
            .def("add_surface", &Scene::add_surface, "shape"_a, "material"_a)
            .def("add_termination", &Scene::add_termination, "shape"_a, "volume"_a = false,
                 "on"_a = Termination::Entry)
            .def("region_at", &Scene::region_at, "point"_a)
            .def("material", &Scene::material, "region"_a)
            .def("starts_terminated", &Scene::starts_terminated, "point"_a)
            .def("next_event", &Scene::next_event, "origin"_a, "direction"_a, "max_distance"_a)
            .def_property_readonly("background", &Scene::background)
            .def_property_readonly("tolerance", &Scene::tolerance)
            .def_property_readonly("objects", &Scene::objects);

    py::class_<Interaction>(m, "Interaction")
            .def_readonly("object", &Interaction::object).def_readonly("kind", &Interaction::kind)
            .def_readonly("before", &Interaction::before).def_readonly("after", &Interaction::after);
    py::class_<State>(m, "State")
            .def(py::init<>())
            .def_readwrite("position", &State::position).def_readwrite("direction", &State::direction)
            .def_readwrite("basis", &State::basis).def_readwrite("curvature", &State::curvature)
            .def_readwrite("field", &State::field).def_readwrite("region", &State::region)
            .def_readwrite("power", &State::power).def_readwrite("distance", &State::distance)
            .def_readwrite("optical_path", &State::opticalPath).def_readonly("history", &State::history)
            .def_readonly("branch", &State::branch).def_readonly("parent", &State::parent);
    m.def("propagate", &propagate, "state"_a, "distance"_a, "material"_a, "frequency_hz"_a);
    py::class_<Segment>(m, "Segment")
            .def_readonly("state", &Segment::state).def_readonly("length", &Segment::length)
            .def_readonly("material", &Segment::material).def_readonly("frequency_hz", &Segment::frequencyHz)
            .def("field_at", &Segment::field_at, "point"_a)
            .def("radii", &Segment::radii, "distance"_a = 0.0);

    py::class_<GaussianBeam>(m, "GaussianBeam")
            .def(py::init([](Vec3 waistPosition, Vec3 direction, const py::object &waistRadius, double power,
                             Vec3C polarization, double launchDistance, double phase,
                             std::optional<Vec3> transverseAxis) {
                     return GaussianBeam{waistPosition, direction, pair(waistRadius), power, polarization,
                                         launchDistance, phase, transverseAxis};
                 }), "waist_position"_a = Vec3::Zero(), "direction"_a = Vec3::UnitZ(), "waist_radius"_a = 0.1,
                 "power"_a = 1.0, "polarization"_a = Vec3C::UnitX(), "launch_distance"_a = 0.0, "phase"_a = 0.0,
                 "transverse_axis"_a = std::nullopt)
            .def_readwrite("waist_position", &GaussianBeam::waistPosition)
            .def_readwrite("direction", &GaussianBeam::direction)
            .def_readwrite("waist_radius", &GaussianBeam::waistRadius)
            .def_readwrite("power", &GaussianBeam::power)
            .def_readwrite("polarization", &GaussianBeam::polarization)
            .def_readwrite("launch_distance", &GaussianBeam::launchDistance)
            .def_readwrite("phase", &GaussianBeam::phase)
            .def_readwrite("transverse_axis", &GaussianBeam::transverseAxis)
            .def("launch", &GaussianBeam::launch, "scene"_a, "frequency_hz"_a)
            .def("beamlets", &GaussianBeam::beamlets, "samples"_a = 17, "radius_ratio"_a = 0.5, "extent"_a = 4.0);
    py::class_<PlaneWave>(m, "PlaneWave")
            .def(py::init<Vec3, Vec3, Vec3C, Complex>(), "position"_a = Vec3::Zero(), "direction"_a = Vec3::UnitZ(),
                 "polarization"_a = Vec3C::UnitX(), "amplitude"_a = 1.0)
            .def_readwrite("position", &PlaneWave::position).def_readwrite("direction", &PlaneWave::direction)
            .def_readwrite("polarization", &PlaneWave::polarization)
            .def_readwrite("amplitude", &PlaneWave::amplitude)
            .def("launch", &PlaneWave::launch, "scene"_a, "frequency_hz"_a);
    py::class_<ApertureSource>(m, "ApertureSource", "Gaussian-beamlet fit of a sampled field on a plane.")
            .def(py::init([](Vec3 position, Vec3 direction, const py::object &halfWidths, int samples,
                             double beamletRadius, std::function<Vec3C(const Vec3 &)> field, double rcond) {
                     return ApertureSource{position, direction, pair(halfWidths), samples, beamletRadius,
                                           std::move(field), rcond, std::nullopt};
                 }), "position"_a, "direction"_a, "half_widths"_a, "samples"_a, "beamlet_radius"_a, "field"_a,
                 "rcond"_a = 1e-10)
            .def_readwrite("position", &ApertureSource::position)
            .def_readwrite("direction", &ApertureSource::direction)
            .def_readwrite("half_widths", &ApertureSource::halfWidths)
            .def_readwrite("samples", &ApertureSource::samples)
            .def_readwrite("beamlet_radius", &ApertureSource::beamletRadius)
            .def_readwrite("field", &ApertureSource::field)
            .def_readwrite("rcond", &ApertureSource::rcond)
            .def_readonly("frequency_hz", &ApertureSource::frequencyHz)
            .def("launch", &ApertureSource::launch, "scene"_a, "frequency_hz"_a)
            .def_static("from_antenna", [](const Antenna &antenna, Vec3 antennaPosition, Vec3 position,
                                           Vec3 direction, const py::object &halfWidths, int samples,
                                           double beamletRadius, double frequencyHz) {
                return ApertureSource::from_antenna(antenna, antennaPosition, position, direction, pair(halfWidths),
                                                    samples, beamletRadius, frequencyHz);
            }, "antenna"_a, "antenna_position"_a, "position"_a, "direction"_a, "half_widths"_a, "samples"_a,
                        "beamlet_radius"_a, "frequency_hz"_a);

    py::class_<SolverConfig>(m, "SolverConfig")
            .def(py::init([](double frequencyHz, int maxInteractions, double maxDistance, double minPowerFraction,
                             std::size_t maxBranches, double maxLossRatio) {
                     return SolverConfig{frequencyHz, maxInteractions, maxDistance, minPowerFraction, maxBranches,
                                         maxLossRatio};
                 }), "frequency_hz"_a = 77e9, "max_interactions"_a = 16, "max_distance"_a = 20.0,
                 "min_power_fraction"_a = 1e-10, "max_branches"_a = 100000, "max_loss_ratio"_a = 0.01)
            .def_readwrite("frequency_hz", &SolverConfig::frequencyHz)
            .def_readwrite("max_interactions", &SolverConfig::maxInteractions)
            .def_readwrite("max_distance", &SolverConfig::maxDistance)
            .def_readwrite("min_power_fraction", &SolverConfig::minPowerFraction)
            .def_readwrite("max_branches", &SolverConfig::maxBranches)
            .def_readwrite("max_loss_ratio", &SolverConfig::maxLossRatio);
    py::class_<Terminal>(m, "Terminal")
            .def_readonly("branch", &Terminal::branch).def_readonly("reason", &Terminal::reason)
            .def_readonly("position", &Terminal::position).def_readonly("power", &Terminal::power)
            .def_readonly("history", &Terminal::history);
    py::class_<Simulation>(m, "Simulation")
            .def_readonly("config", &Simulation::config).def_readonly("scene", &Simulation::scene)
            .def_readonly("segments", &Simulation::segments).def_readonly("terminals", &Simulation::terminals)
            .def_readonly("warnings", &Simulation::warnings)
            .def_readonly("launched_power", &Simulation::launchedPower)
            .def_readonly("absorbed_power", &Simulation::absorbedPower)
            .def_readonly("interface_flux_residual", &Simulation::interfaceFluxResidual)
            .def("field_at", &Simulation::field_at, "point"_a, py::call_guard<py::gil_scoped_release>())
            .def("impulse_response", &Simulation::impulse_response, "point"_a, "analyzer"_a = Vec3C::UnitX(),
                 "delay_tolerance"_a = 1e-13, py::call_guard<py::gil_scoped_release>())
            .def("fields_at", [](const Simulation &s, const Points &points) {
                py::gil_scoped_release release;
                Fields values(points.rows(), 3);
                for (Eigen::Index i = 0; i < points.rows(); ++i)
                    values.row(i) = s.field_at(points.row(i).transpose()).transpose();
                return values;
            }, "points"_a, "Complex E at each row of an (N, 3) array.")
            .def("power_balance", &Simulation::power_balance);

    // Sources run under the GIL (aperture callbacks may be Python); transport runs without it.
    m.def("solve", [](const Scene &scene, const py::object &source, const SolverConfig &config) {
        std::vector<State> states;
        for (const auto &s: sources(source)) {
            auto launched = launch(scene, s, config.frequencyHz);
            states.insert(states.end(), launched.begin(), launched.end());
        }
        py::gil_scoped_release release;
        return solve_states(scene, std::move(states), config);
    }, "scene"_a, "sources"_a, "config"_a = SolverConfig{}, "Launch one source or a list of sources and trace them.");
    m.def("solve_states", &solve_states, "scene"_a, "states"_a, "config"_a = SolverConfig{},
          py::call_guard<py::gil_scoped_release>());
}
