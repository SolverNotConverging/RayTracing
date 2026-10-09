#include "Solver.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rt {
    namespace {
        bool positive(double x) { return std::isfinite(x) && x > 0; }

        void validate_config(const SolverConfig &s) {
            if (!positive(s.frequencyHz) || s.rayCount == 0 || s.maxReflections < 0 ||
                !positive(s.maxDistance) || !positive(s.receptionRadius) ||
                !positive(s.commonSourceReference) ||
                !std::isfinite(s.delayToleranceSeconds) || s.delayToleranceSeconds < 0)
                throw std::invalid_argument(
                    "Invalid solver frequency, ray count, limits or normalization");
            s.medium.refractive_index();
            // Exercise the low-level option validators before any rays are launched.
            auto refine = s.refinement;
            refine.maxPathDistance_ = s.maxDistance;
            refine_path(Vec3::Zero(), Vec3::UnitX(), Vec3::UnitX(), {}, {}, refine);
            calculate_spreading(Vec3::Zero(), Vec3::UnitX(), Vec3::UnitX(), {}, {},
                                s.spreading);
            deduplicate_paths({}, s.deduplication);
        }

        void make_response(SimulationResult &result) {
            const auto &cfg = result.settings;
            result.receiver.antenna.validate(cfg.frequencyHz, cfg.medium);
            // Resolve receive normalization now, including missing imported input power.
            result.receiver.antenna.reference_power(cfg.frequencyHz, cfg.medium);
            std::vector<ReconstructedField> projected;
            result.responseRayIndices.clear();
            for (std::size_t i = 0; i < result.rays.size(); ++i) {
                auto &ray = result.rays[i];
                if (!ray.field)
                    continue;
                ray.receivedCoefficient = result.receiver.antenna.receive(
                    -ray.field->arrivalDirection_, ray.field->normalizedReceiverField_,
                    cfg.frequencyHz, cfg.medium);
                ReconstructedField scalar = *ray.field;
                scalar.normalizedReceiverField_ =
                        ray.receivedCoefficient * Vec3::UnitZ().cast<Complex>();
                projected.push_back(std::move(scalar));
                result.responseRayIndices.push_back(i);
            }
            result.impulseResponse = calculate_impulse_response(
                projected, Vec3::UnitZ().cast<Complex>(), cfg.delayToleranceSeconds);
            result.impulseResponse.frequencyHz_ = cfg.frequencyHz;
            result.impulseResponse.receiverPolarization_ = Vec3C::Zero();
            result.impulseResponse.receiverModel_ =
                    "reciprocal antenna port (kind " +
                    std::to_string(static_cast<int>(result.receiver.antenna.kind)) + ")";
            // Restore actual incident vector sums: scalar projection and vector data are
            // separate.
            for (auto &tap: result.impulseResponse.taps_) {
                tap.normalizedField_ = Vec3C::Zero();
                for (auto index: tap.pathIndices_)
                    tap.normalizedField_ +=
                            result.rays.at(result.responseRayIndices.at(index))
                            .field->normalizedReceiverField_;
            }
        }
    } // namespace

    SimulationResult solve(const Scene &scene, const Transmitter &tx,
                           const Receiver &rx, const SolverConfig &settings) {
        validate_config(settings);
        if (!tx.position.allFinite() || !rx.position.allFinite() ||
            (tx.position - rx.position).norm() == 0)
            throw std::invalid_argument("Tx and Rx must be finite and distinct");
        tx.antenna.validate(settings.frequencyHz, settings.medium);
        rx.antenna.validate(settings.frequencyHz, settings.medium);
        rx.antenna.reference_power(settings.frequencyHz, settings.medium);
        SimulationResult result;
        result.settings = settings;
        result.scene = scene;
        result.transmitter = tx;
        result.receiver = rx;
        TraceOptions tracing{
            settings.maxReflections, settings.maxDistance,
            settings.medium.refractive_index()
        };
        auto refinement = settings.refinement;
        refinement.maxPathDistance_ = settings.maxDistance;
        auto spread = settings.spreading;
        spread.receiverTolerance_ = refinement.receiverTolerance_;
        spread.tMin_ = refinement.tMin_;
        // Include the exact direct direction as well as sampled rays; deduplication
        // removes repeats, and obstacles still participate in every trace.
        auto directions = launch_directions(settings.rayCount);
        directions.push_back((rx.position - tx.position).normalized());
        std::vector<RefinementResult> trials;
        for (const auto &direction: directions) {
            // Geometry discovery is independent of source-pattern nulls.
            EMRay ray(tx.position, direction, launch_polarization(direction),
                      C / settings.frequencyHz);
            const auto traced = trace_ray(ray, scene.surfaces(), rx.position,
                                          settings.receptionRadius, tracing);
            ++result.launchedRays;
            ++result.traceStatusCounts.at(static_cast<std::size_t>(traced.status_));
            if (traced.status_ != TraceStatus::Received)
                continue;
            std::vector<std::size_t> sequence;
            for (const auto &hit: ray.reflections_)
                sequence.push_back(hit.surfaceIndex_);
            auto refined = refine_path(tx.position, direction, rx.position,
                                       scene.surfaces(), sequence, refinement);
            Candidate candidate;
            candidate.refinement = refined;
            for (const auto &point: ray.path_)
                candidate.coarseVertices.push_back(point.position_);
            result.candidates.push_back(std::move(candidate));
            trials.push_back(std::move(refined));
        }
        for (auto index: deduplicate_paths(trials, settings.deduplication)) {
            const auto &path = trials[index];
            std::vector<std::size_t> sequence;
            for (const auto &hit: path.geometry_.reflections_)
                sequence.push_back(hit.surfaceIndex_);
            SolvedRay ray;
            ray.candidateIndex = index;
            ray.refinement = path;
            ray.spreading =
                    calculate_spreading(tx.position, path.launchDirection_, rx.position,
                                        scene.surfaces(), sequence, spread);
            ray.sourceField =
                    tx.antenna.farfield(path.launchDirection_, settings.frequencyHz,
                                        settings.medium) /
                    spread.referenceDistance_;
            if (ray.spreading.status_ == SpreadingStatus::Valid) {
                FieldReconstructionOptions field;
                field.frequencyHz_ = settings.frequencyHz;
                field.refractiveIndex_ = tracing.refractiveIndex_;
                field.receiverTolerance_ = refinement.receiverTolerance_;
                field.sourceReferenceAmplitude_ = settings.commonSourceReference;
                ray.field =
                        reconstruct_field(path, ray.spreading, ray.sourceField, field);
            }
            result.rays.push_back(std::move(ray));
        }
        make_response(result);
        return result;
    }

    void update_receiver(SimulationResult &result, const Antenna &antenna) {
        // Validate and calculate in a copy so a failed antenna update is reversible.
        SimulationResult changed = result;
        changed.receiver.antenna = antenna;
        make_response(changed);
        result = std::move(changed);
    }
} // namespace rt
