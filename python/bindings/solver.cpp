#include "bindings.hpp"
#include <RayTracing/ResultIO.hpp>

void bind_solver(py::module_ &m) {
    using namespace rt;
    py::class_<RefinementOptions>(m, "RefinementOptions")
            .def(py::init<>()).def_readwrite("max_iterations", &RefinementOptions::maxIterations_)
            .def_readwrite("receiver_tolerance", &RefinementOptions::receiverTolerance_)
            .def_readwrite("finite_difference_step", &RefinementOptions::finiteDifferenceStep_)
            .def_readwrite("max_direction_step", &RefinementOptions::maxDirectionStep_)
            .def_readwrite("max_backtracks", &RefinementOptions::maxBacktracks_)
            .def_readwrite("t_min", &RefinementOptions::tMin_)
            .def_readwrite("max_path_distance", &RefinementOptions::maxPathDistance_)
            .def_readwrite("corner_separation_tolerance", &RefinementOptions::cornerSeparationTolerance_);
    py::class_<SpreadingOptions>(m, "SpreadingOptions")
            .def(py::init<>()).def_readwrite("angular_step", &SpreadingOptions::angularStep_)
            .def_readwrite("max_step_halvings", &SpreadingOptions::maxStepHalvings_)
            .def_readwrite("receiver_tolerance", &SpreadingOptions::receiverTolerance_)
            .def_readwrite("relative_derivative_tolerance", &SpreadingOptions::relativeDerivativeTolerance_)
            .def_readwrite("absolute_derivative_tolerance", &SpreadingOptions::absoluteDerivativeTolerance_)
            .def_readwrite("minimum_singular_value", &SpreadingOptions::minimumSingularValue_)
            .def_readwrite("minimum_singular_value_ratio", &SpreadingOptions::minimumSingularValueRatio_)
            .def_readwrite("reference_distance", &SpreadingOptions::referenceDistance_)
            .def_readwrite("t_min", &SpreadingOptions::tMin_);
    py::class_<DeduplicationOptions>(m, "DeduplicationOptions")
            .def(py::init<>()).def_readwrite("position_tolerance", &DeduplicationOptions::positionTolerance_)
            .def_readwrite("angle_tolerance", &DeduplicationOptions::angleTolerance_)
            .def_readwrite("length_tolerance", &DeduplicationOptions::lengthTolerance_);
    py::class_<SolverConfig>(m, "SolverConfig")
            .def(py::init([](double frequency, std::size_t count, int reflections, double distance, double radius) {
                     SolverConfig c;
                     c.frequencyHz = frequency;
                     c.rayCount = count;
                     c.maxReflections = reflections;
                     c.maxDistance = distance;
                     c.receptionRadius = radius;
                     return c;
                 }), py::arg("frequency_hz") = 77e9, py::arg("ray_count") = 2000,
                 py::arg("max_reflections") = 8, py::arg("max_distance") = 20.0, py::arg("reception_radius") = 0.12)
            .def_readwrite("frequency_hz", &SolverConfig::frequencyHz).def_readwrite(
                "ray_count", &SolverConfig::rayCount)
            .def_readwrite("max_reflections", &SolverConfig::maxReflections).def_readwrite(
                "max_distance", &SolverConfig::maxDistance)
            .def_readwrite("reception_radius", &SolverConfig::receptionRadius)
            .def_readwrite("common_source_reference", &SolverConfig::commonSourceReference)
            .def_readwrite("delay_tolerance_seconds", &SolverConfig::delayToleranceSeconds)
            .def_readwrite("medium", &SolverConfig::medium).def_readwrite("refinement", &SolverConfig::refinement)
            .def_readwrite("spreading", &SolverConfig::spreading).def_readwrite(
                "deduplication", &SolverConfig::deduplication);

    // Snapshot mutable inputs under the GIL, then let Python remain responsive while solving.
    m.def("solve", [](Scene scene, Transmitter tx, Receiver rx, SolverConfig config) {
        py::gil_scoped_release release;
        return solve(scene, tx, rx, config);
    }, py::arg("scene"), py::arg("transmitter"), py::arg("receiver"), py::arg("config") = SolverConfig{});
    m.def("update_receiver", &update_receiver, py::arg("result"), py::arg("antenna"));
    m.def("save_h5", &save_h5, py::arg("result"), py::arg("path"), py::call_guard<py::gil_scoped_release>());
    m.def("load_h5", &load_h5, py::arg("path"), py::call_guard<py::gil_scoped_release>());
    m.def("save_impulse_csv", &save_impulse_csv, py::arg("response"), py::arg("path"));
    m.def("load_impulse_csv", &load_impulse_csv, py::arg("path"));
    m.def("frequency_response", &frequency_response, py::arg("response"), py::arg("offset_hz") = 0.0);
    m.def("evaluate_sequence", &evaluate_sequence, py::arg("transmitter_position"), py::arg("launch_direction"),
          py::arg("receiver_position"), py::arg("surfaces"), py::arg("surface_sequence"), py::arg("t_min") = 1e-8);
    m.def("refine_path", &refine_path, py::arg("transmitter_position"), py::arg("launch_direction"),
          py::arg("receiver_position"), py::arg("surfaces"), py::arg("surface_sequence"),
          py::arg("options") = RefinementOptions{});
    m.def("calculate_spreading", &calculate_spreading, py::arg("transmitter_position"), py::arg("launch_direction"),
          py::arg("receiver_position"), py::arg("surfaces"), py::arg("surface_sequence"),
          py::arg("options") = SpreadingOptions{});
    m.def("deduplicate_paths", &deduplicate_paths, py::arg("paths"), py::arg("options") = DeduplicationOptions{});
}
