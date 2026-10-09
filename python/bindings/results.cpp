#include "bindings.hpp"

void bind_results(py::module_ &m) {
    py::enum_<SequenceStatus>(m, "SequenceStatus")
            .value("VALID", SequenceStatus::Valid)
            .value("MISSED_SURFACE", SequenceStatus::MissedSurface)
            .value("UNEXPECTED_SURFACE", SequenceStatus::UnexpectedSurface)
            .value("RECEIVER_BEHIND", SequenceStatus::ReceiverBehind)
            .value("FINAL_SEGMENT_BLOCKED", SequenceStatus::FinalSegmentBlocked)
            .value("AMBIGUOUS_HIT", SequenceStatus::AmbiguousHit);
    py::enum_<RefinementStatus>(m, "RefinementStatus")
            .value("CONVERGED", RefinementStatus::Converged)
            .value("INVALID_GEOMETRY", RefinementStatus::InvalidGeometry)
            .value("DERIVATIVE_UNAVAILABLE", RefinementStatus::DerivativeUnavailable)
            .value("SINGULAR_JACOBIAN", RefinementStatus::SingularJacobian)
            .value("STALLED", RefinementStatus::Stalled)
            .value("ITERATION_LIMIT", RefinementStatus::IterationLimit)
            .value("PATH_TOO_LONG", RefinementStatus::PathTooLong)
            .value("UNRESOLVED_CORNER", RefinementStatus::UnresolvedCorner)
            .value("UNRESOLVED_EDGE", RefinementStatus::UnresolvedEdge);
    py::enum_<SpreadingStatus>(m, "SpreadingStatus")
            .value("VALID", SpreadingStatus::Valid)
            .value("INVALID_CENTRAL_PATH", SpreadingStatus::InvalidCentralPath)
            .value("INVALID_PERTURBATION", SpreadingStatus::InvalidPerturbation)
            .value("UNSTABLE_DERIVATIVE", SpreadingStatus::UnstableDerivative)
            .value("CAUSTIC", SpreadingStatus::Caustic);
    py::enum_<TraceStatus>(m, "TraceStatus")
            .value("RECEIVED", TraceStatus::Received)
            .value("ESCAPED", TraceStatus::Escaped)
            .value("REFLECTION_LIMIT", TraceStatus::ReflectionLimit)
            .value("DISTANCE_LIMIT", TraceStatus::DistanceLimit)
            .value("AMBIGUOUS_HIT", TraceStatus::AmbiguousHit);
    py::class_<SequenceReflection>(m, "SequenceReflection")
            .def_readonly("surface_index", &SequenceReflection::surfaceIndex_)
            .def_readonly("position", &SequenceReflection::position_)
            .def_readonly("normal", &SequenceReflection::normal_)
            .def_readonly("segment_distance", &SequenceReflection::segmentDistance_);
    py::class_<ReceiverApproach>(m, "ReceiverApproach")
            .def_readonly("closest_point", &ReceiverApproach::closestPoint_)
            .def_readonly("residual", &ReceiverApproach::residual_)
            .def_readonly("miss_distance", &ReceiverApproach::missDistance_)
            .def_readonly("final_segment_distance", &ReceiverApproach::finalSegmentDistance_)
            .def_readonly("path_distance", &ReceiverApproach::pathDistance_);
    py::class_<SequenceEvaluation>(m, "SequenceEvaluation")
            .def_readonly("status", &SequenceEvaluation::status_)
            .def_readonly("reflections", &SequenceEvaluation::reflections_)
            .def_readonly("final_origin", &SequenceEvaluation::finalOrigin_)
            .def_readonly("final_direction", &SequenceEvaluation::finalDirection_)
            .def_readonly("receiver", &SequenceEvaluation::receiver_)
            .def_readonly("blocking_surface_index", &SequenceEvaluation::blockingSurfaceIndex_);
    py::class_<RefinementResult>(m, "RefinementResult")
            .def_readonly("status", &RefinementResult::status_)
            .def_readonly("transmitter_position", &RefinementResult::transmitterPosition_)
            .def_readonly("receiver_position", &RefinementResult::receiverPosition_)
            .def_readonly("launch_direction", &RefinementResult::launchDirection_)
            .def_readonly("geometry", &RefinementResult::geometry_)
            .def_readonly("iterations", &RefinementResult::iterations_);
    py::class_<SpreadingResult>(m, "SpreadingResult")
            .def_readonly("status", &SpreadingResult::status_)
            .def_readonly("jacobian", &SpreadingResult::jacobian_)
            .def_readonly("singular_values", &SpreadingResult::singularValues_)
            .def_readonly("area_per_solid_angle", &SpreadingResult::areaPerSolidAngle_)
            .def_readonly("angular_step", &SpreadingResult::angularStep_)
            .def_readonly("derivative_difference", &SpreadingResult::derivativeDifference_)
            .def_readonly("field_factor", &SpreadingResult::fieldFactor_)
            .def_readonly("launch_u", &SpreadingResult::launchU_)
            .def_readonly("launch_v", &SpreadingResult::launchV_)
            .def_readonly("receiver_u", &SpreadingResult::receiverU_)
            .def_readonly("receiver_v", &SpreadingResult::receiverV_);
    py::class_<ReconstructedField>(m, "ReconstructedField")
            .def_readonly("receiver_field", &ReconstructedField::receiverField_)
            .def_readonly("normalized_receiver_field", &ReconstructedField::normalizedReceiverField_)
            .def_readonly("transported_reference_field", &ReconstructedField::transportedReferenceField_)
            .def_readonly("arrival_direction", &ReconstructedField::arrivalDirection_)
            .def_readonly("path_distance", &ReconstructedField::pathDistance_)
            .def_readonly("optical_path", &ReconstructedField::opticalPath_)
            .def_readonly("delay_seconds", &ReconstructedField::delaySeconds_)
            .def_readonly("frequency_hz", &ReconstructedField::frequencyHz_)
            .def_readonly("field_factor", &ReconstructedField::fieldFactor_);
    py::class_<ImpulseTap>(m, "ImpulseTap")
            .def_readonly("delay_seconds", &ImpulseTap::delaySeconds_)
            .def_readonly("coefficient", &ImpulseTap::coefficient_)
            .def_readonly("normalized_field", &ImpulseTap::normalizedField_)
            .def_readonly("path_indices", &ImpulseTap::pathIndices_);
    py::class_<ImpulseResponse>(m, "ImpulseResponse")
            .def_readonly("frequency_hz", &ImpulseResponse::frequencyHz_)
            .def_readonly("receiver_polarization", &ImpulseResponse::receiverPolarization_)
            .def_readonly("taps", &ImpulseResponse::taps_)
            .def_readonly("delay_tolerance_seconds", &ImpulseResponse::delayToleranceSeconds_)
            .def_readonly("receiver_model", &ImpulseResponse::receiverModel_);
    py::class_<rt::Candidate>(m, "Candidate")
            .def_readonly("refinement", &rt::Candidate::refinement)
            .def_readonly("coarse_vertices", &rt::Candidate::coarseVertices);
    py::class_<rt::SolvedRay>(m, "SolvedRay")
            .def_readonly("candidate_index", &rt::SolvedRay::candidateIndex)
            .def_readonly("refinement", &rt::SolvedRay::refinement)
            .def_readonly("spreading", &rt::SolvedRay::spreading)
            .def_readonly("field", &rt::SolvedRay::field)
            .def_readonly("source_field", &rt::SolvedRay::sourceField)
            .def_readonly("received_coefficient", &rt::SolvedRay::receivedCoefficient)
            .def_property_readonly("vertices", [](const rt::SolvedRay &ray) {
                const auto &path = ray.refinement;
                const auto &geometry = path.geometry_;
                const auto count = 1 + geometry.reflections_.size() + (geometry.receiver_ ? 1 : 0);
                Eigen::Matrix<double, Eigen::Dynamic, 3> points(count, 3);
                points.row(0) = path.transmitterPosition_.transpose();
                for (std::size_t i = 0; i < geometry.reflections_.size(); ++i)
                    points.row(i + 1) = geometry.reflections_[i].position_.transpose();
                if (geometry.receiver_) points.row(count - 1) = geometry.receiver_->closestPoint_.transpose();
                return points;
            });
    py::class_<rt::SimulationResult>(m, "SimulationResult")
            .def_readonly("settings", &rt::SimulationResult::settings)
            .def_readonly("scene", &rt::SimulationResult::scene)
            .def_readonly("transmitter", &rt::SimulationResult::transmitter)
            .def_readonly("receiver", &rt::SimulationResult::receiver)
            .def_readonly("rays", &rt::SimulationResult::rays)
            .def_readonly("candidates", &rt::SimulationResult::candidates)
            .def_readonly("impulse_response", &rt::SimulationResult::impulseResponse)
            .def_readonly("launched_rays", &rt::SimulationResult::launchedRays)
            .def_readonly("trace_status_counts", &rt::SimulationResult::traceStatusCounts)
            .def_readonly("response_ray_indices", &rt::SimulationResult::responseRayIndices);
}
