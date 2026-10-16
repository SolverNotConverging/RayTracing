// Point-source multipath solver.
#include "bindings.hpp"

using namespace rt;

namespace {
    using Points = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

    Points to_points(const std::vector<Vec3> &vertices) {
        Points points(vertices.size(), 3);
        for (std::size_t i = 0; i < vertices.size(); ++i) points.row(i) = vertices[i];
        return points;
    }
} // namespace

void bind_point_ray(py::module_ &m) {
    py::class_<Scene>(m, "Scene", "Reflecting surfaces for the point-ray solver.")
            .def(py::init<>())
            .def("add", &Scene::add, "surface"_a, "material"_a = Material::pec(),
                 "Add a surface; returns its index in paths.")
            .def_property_readonly("surfaces", &Scene::surfaces)
            .def_property_readonly("materials", &Scene::materials)
            .def("__len__", [](const Scene &s) { return s.surfaces().size(); });
    py::class_<Transmitter>(m, "Transmitter")
            .def(py::init<Vec3, Antenna>(), "position"_a, "antenna"_a = Antenna{})
            .def_readwrite("position", &Transmitter::position).def_readwrite("antenna", &Transmitter::antenna);
    py::class_<Receiver>(m, "Receiver")
            .def(py::init<Vec3, Antenna>(), "position"_a, "antenna"_a = Antenna{})
            .def_readwrite("position", &Receiver::position).def_readwrite("antenna", &Receiver::antenna);

    py::class_<RefinementOptions>(m, "RefinementOptions")
            .def(py::init<>())
            .def_readwrite("max_iterations", &RefinementOptions::maxIterations)
            .def_readwrite("receiver_tolerance", &RefinementOptions::receiverTolerance)
            .def_readwrite("finite_difference_step", &RefinementOptions::finiteDifferenceStep)
            .def_readwrite("max_direction_step", &RefinementOptions::maxDirectionStep)
            .def_readwrite("max_backtracks", &RefinementOptions::maxBacktracks)
            .def_readwrite("t_min", &RefinementOptions::tMin)
            .def_readwrite("corner_separation_tolerance", &RefinementOptions::cornerSeparationTolerance);
    py::class_<SpreadingOptions>(m, "SpreadingOptions")
            .def(py::init<>())
            .def_readwrite("angular_step", &SpreadingOptions::angularStep)
            .def_readwrite("max_step_halvings", &SpreadingOptions::maxStepHalvings)
            .def_readwrite("relative_derivative_tolerance", &SpreadingOptions::relativeDerivativeTolerance)
            .def_readwrite("absolute_derivative_tolerance", &SpreadingOptions::absoluteDerivativeTolerance)
            .def_readwrite("minimum_singular_value", &SpreadingOptions::minimumSingularValue)
            .def_readwrite("minimum_singular_value_ratio", &SpreadingOptions::minimumSingularValueRatio)
            .def_readwrite("reference_distance", &SpreadingOptions::referenceDistance);
    py::class_<DeduplicationOptions>(m, "DeduplicationOptions")
            .def(py::init<>())
            .def_readwrite("position_tolerance", &DeduplicationOptions::positionTolerance)
            .def_readwrite("angle_tolerance", &DeduplicationOptions::angleTolerance)
            .def_readwrite("length_tolerance", &DeduplicationOptions::lengthTolerance);
    py::class_<SolverConfig>(m, "SolverConfig")
            .def(py::init([](double frequencyHz, std::size_t rayCount, int maxReflections, double maxDistance,
                             double receptionRadius, Material background) {
                     SolverConfig c;
                     c.frequencyHz = frequencyHz;
                     c.rayCount = rayCount;
                     c.maxReflections = maxReflections;
                     c.maxDistance = maxDistance;
                     c.receptionRadius = receptionRadius;
                     c.background = std::move(background);
                     return c;
                 }), "frequency_hz"_a = 77e9, "ray_count"_a = 2000, "max_reflections"_a = 8, "max_distance"_a = 20.0,
                 "reception_radius"_a = 0.12, "background"_a = Material{})
            .def_readwrite("frequency_hz", &SolverConfig::frequencyHz)
            .def_readwrite("ray_count", &SolverConfig::rayCount)
            .def_readwrite("max_reflections", &SolverConfig::maxReflections)
            .def_readwrite("max_distance", &SolverConfig::maxDistance)
            .def_readwrite("reception_radius", &SolverConfig::receptionRadius)
            .def_readwrite("source_reference", &SolverConfig::sourceReference)
            .def_readwrite("delay_tolerance", &SolverConfig::delayTolerance)
            .def_readwrite("background", &SolverConfig::background)
            .def_readwrite("refinement", &SolverConfig::refinement)
            .def_readwrite("spreading", &SolverConfig::spreading)
            .def_readwrite("deduplication", &SolverConfig::deduplication);

    py::enum_<TraceStatus>(m, "TraceStatus")
            .value("RECEIVED", TraceStatus::Received).value("ESCAPED", TraceStatus::Escaped)
            .value("REFLECTION_LIMIT", TraceStatus::ReflectionLimit)
            .value("DISTANCE_LIMIT", TraceStatus::DistanceLimit).value("AMBIGUOUS_HIT", TraceStatus::AmbiguousHit);
    py::enum_<SequenceStatus>(m, "SequenceStatus")
            .value("VALID", SequenceStatus::Valid).value("MISSED_SURFACE", SequenceStatus::MissedSurface)
            .value("UNEXPECTED_SURFACE", SequenceStatus::UnexpectedSurface)
            .value("RECEIVER_BEHIND", SequenceStatus::ReceiverBehind)
            .value("FINAL_SEGMENT_BLOCKED", SequenceStatus::FinalSegmentBlocked)
            .value("AMBIGUOUS_HIT", SequenceStatus::AmbiguousHit);
    py::enum_<RefinementStatus>(m, "RefinementStatus")
            .value("CONVERGED", RefinementStatus::Converged)
            .value("INVALID_GEOMETRY", RefinementStatus::InvalidGeometry)
            .value("DERIVATIVE_UNAVAILABLE", RefinementStatus::DerivativeUnavailable)
            .value("SINGULAR_JACOBIAN", RefinementStatus::SingularJacobian)
            .value("STALLED", RefinementStatus::Stalled).value("ITERATION_LIMIT", RefinementStatus::IterationLimit)
            .value("PATH_TOO_LONG", RefinementStatus::PathTooLong)
            .value("UNRESOLVED_CORNER", RefinementStatus::UnresolvedCorner)
            .value("UNRESOLVED_EDGE", RefinementStatus::UnresolvedEdge);
    py::enum_<SpreadingStatus>(m, "SpreadingStatus")
            .value("VALID", SpreadingStatus::Valid).value("INVALID_CENTRAL_PATH", SpreadingStatus::InvalidCentralPath)
            .value("INVALID_PERTURBATION", SpreadingStatus::InvalidPerturbation)
            .value("UNSTABLE_DERIVATIVE", SpreadingStatus::UnstableDerivative)
            .value("CAUSTIC", SpreadingStatus::Caustic);

    py::class_<SequenceReflection>(m, "SequenceReflection")
            .def_readonly("surface_index", &SequenceReflection::surfaceIndex)
            .def_readonly("position", &SequenceReflection::position)
            .def_readonly("normal", &SequenceReflection::normal)
            .def_readonly("segment_distance", &SequenceReflection::segmentDistance);
    py::class_<ReceiverApproach>(m, "ReceiverApproach")
            .def_readonly("closest_point", &ReceiverApproach::closestPoint)
            .def_readonly("miss_distance", &ReceiverApproach::missDistance)
            .def_readonly("path_distance", &ReceiverApproach::pathDistance);
    py::class_<SequenceEvaluation>(m, "SequenceEvaluation")
            .def_readonly("status", &SequenceEvaluation::status)
            .def_readonly("reflections", &SequenceEvaluation::reflections)
            .def_readonly("final_direction", &SequenceEvaluation::finalDirection)
            .def_readonly("receiver", &SequenceEvaluation::receiver);
    py::class_<RefinementResult>(m, "RefinementResult")
            .def_readonly("status", &RefinementResult::status)
            .def_readonly("launch_direction", &RefinementResult::launchDirection)
            .def_readonly("geometry", &RefinementResult::geometry)
            .def_readonly("iterations", &RefinementResult::iterations);
    py::class_<SpreadingResult>(m, "SpreadingResult")
            .def_readonly("status", &SpreadingResult::status)
            .def_readonly("jacobian", &SpreadingResult::jacobian)
            .def_readonly("singular_values", &SpreadingResult::singularValues)
            .def_readonly("field_factor", &SpreadingResult::fieldFactor)
            .def_readonly("caustic_count", &SpreadingResult::causticCount);
    py::class_<Reflection>(m, "Reflection")
            .def_readonly("incidence_angle", &Reflection::incidenceAngle)
            .def_readonly("s", &Reflection::s).def_readonly("p", &Reflection::p);
    py::class_<PathField>(m, "PathField")
            .def_readonly("field", &PathField::field)
            .def_readonly("normalized_field", &PathField::normalizedField)
            .def_readonly("arrival_direction", &PathField::arrivalDirection)
            .def_readonly("reflections", &PathField::reflections)
            .def_readonly("path_length", &PathField::pathLength)
            .def_readonly("delay", &PathField::delay)
            .def_readonly("caustic_count", &PathField::causticCount);
    py::class_<ImpulseTap>(m, "ImpulseTap")
            .def_readonly("delay", &ImpulseTap::delay)
            .def_readonly("coefficient", &ImpulseTap::coefficient)
            .def_readonly("field", &ImpulseTap::field)
            .def_readonly("paths", &ImpulseTap::paths);
    py::class_<ImpulseResponse>(m, "ImpulseResponse")
            .def_readonly("frequency_hz", &ImpulseResponse::frequencyHz)
            .def_readonly("taps", &ImpulseResponse::taps)
            .def_property_readonly("delays", [](const ImpulseResponse &r) {
                Eigen::VectorXd values(r.taps.size());
                for (std::size_t i = 0; i < r.taps.size(); ++i) values[i] = r.taps[i].delay;
                return values;
            })
            .def_property_readonly("coefficients", [](const ImpulseResponse &r) {
                Eigen::VectorXcd values(r.taps.size());
                for (std::size_t i = 0; i < r.taps.size(); ++i) values[i] = r.taps[i].coefficient;
                return values;
            })
            .def("frequency_response", &frequency_response, "offset_hz"_a = 0.0);
    py::class_<Candidate>(m, "Candidate")
            .def_readonly("refinement", &Candidate::refinement)
            .def_property_readonly("coarse_vertices", [](const Candidate &c) { return to_points(c.coarseVertices); });
    py::class_<Path>(m, "Path")
            .def_readonly("candidate", &Path::candidate)
            .def_readonly("refinement", &Path::refinement)
            .def_readonly("spreading", &Path::spreading)
            .def_readonly("field", &Path::field)
            .def_readonly("source_field", &Path::sourceField)
            .def_readonly("coefficient", &Path::coefficient)
            .def_property_readonly("vertices", [](const Path &p) { return to_points(p.vertices()); });
    py::class_<SimulationResult>(m, "SimulationResult")
            .def_readonly("config", &SimulationResult::config)
            .def_readonly("scene", &SimulationResult::scene)
            .def_readonly("transmitter", &SimulationResult::transmitter)
            .def_readonly("receiver", &SimulationResult::receiver)
            .def_readonly("paths", &SimulationResult::paths)
            .def_readonly("candidates", &SimulationResult::candidates)
            .def_readonly("impulse_response", &SimulationResult::impulseResponse)
            .def_readonly("trace_status_counts", &SimulationResult::traceStatusCounts);

    m.def("solve", &solve, "scene"_a, "transmitter"_a, "receiver"_a, "config"_a = SolverConfig{},
          py::call_guard<py::gil_scoped_release>());
    m.def("set_receiver_antenna", &set_receiver_antenna, "result"_a, "antenna"_a,
          "Reapply a receive antenna to a result without retracing.");
    m.def("launch_directions", &launch_directions, "count"_a);
    m.def("refine_path", &refine_path, "transmitter_position"_a, "launch_direction"_a, "receiver_position"_a,
          "surfaces"_a, "surface_sequence"_a, "options"_a = RefinementOptions{});
    m.def("calculate_spreading", &calculate_spreading, "transmitter_position"_a, "launch_direction"_a,
          "receiver_position"_a, "surfaces"_a, "surface_sequence"_a, "options"_a = SpreadingOptions{});
}
