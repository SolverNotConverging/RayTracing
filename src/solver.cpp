#include <raytracing/solver.hpp>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace rt {
    namespace {
        void require(bool condition, const char *message) {
            if (!condition) throw std::invalid_argument(message);
        }

        bool positive(double x) { return std::isfinite(x) && x > 0; }

        void validate(const SolverConfig &c) {
            require(positive(c.frequencyHz) && c.rayCount > 0 && c.maxReflections >= 0 && positive(c.maxDistance) &&
                    positive(c.receptionRadius) && positive(c.sourceReference) && c.delayTolerance >= 0,
                    "Invalid solver frequency, ray count, limits or normalization");
            c.background.validate();
            require(c.background.kind == MaterialKind::Dielectric && c.background.index(c.frequencyHz).imag() == 0,
                    "The point-ray background must be a lossless dielectric");
            require(positive(c.refinement.receiverTolerance) && positive(c.refinement.finiteDifferenceStep) &&
                    c.spreading.angularStep > 0 && c.spreading.angularStep <= 0.1 &&
                    positive(c.spreading.referenceDistance) && positive(c.deduplication.positionTolerance) &&
                    positive(c.deduplication.angleTolerance) && positive(c.deduplication.lengthTolerance),
                    "Invalid refinement, spreading or deduplication options");
        }

        bool curved(const Surface &surface) {
            return std::holds_alternative<Sphere>(surface) || std::holds_alternative<Cylinder>(surface);
        }

        std::vector<std::size_t> sequence_of(const SequenceEvaluation &geometry) {
            std::vector<std::size_t> sequence;
            for (const auto &hit: geometry.reflections) sequence.push_back(hit.surfaceIndex);
            return sequence;
        }

        void make_response(SimulationResult &result) {
            const auto &c = result.config;
            const auto &antenna = result.receiver.antenna;
            antenna.validate(c.frequencyHz, c.background);
            std::vector<ImpulseTap> paths;
            for (std::size_t i = 0; i < result.paths.size(); ++i) {
                auto &path = result.paths[i];
                if (!path.field) continue;
                path.coefficient = antenna.receive(-path.field->arrivalDirection, path.field->normalizedField,
                                                   c.frequencyHz, c.background);
                paths.push_back({path.field->delay, path.coefficient, path.field->normalizedField, {i}});
            }
            result.impulseResponse = group_taps(std::move(paths), c.frequencyHz, c.delayTolerance);
        }
    } // namespace

    std::size_t Scene::add(Surface surface, Material material) {
        validate_surface(surface);
        material.validate();
        surfaces_.push_back(std::move(surface));
        materials_.push_back(std::move(material));
        return surfaces_.size() - 1;
    }

    std::vector<Vec3> Path::vertices() const {
        const auto &geometry = refinement.geometry;
        std::vector<Vec3> points{refinement.transmitterPosition};
        for (const auto &hit: geometry.reflections) points.push_back(hit.position);
        if (geometry.receiver) points.push_back(geometry.receiver->closestPoint);
        return points;
    }

    SimulationResult solve(const Scene &scene, const Transmitter &tx, const Receiver &rx, const SolverConfig &config) {
        validate(config);
        require(tx.position.allFinite() && rx.position.allFinite(), "Tx and Rx positions must be finite");
        tx.antenna.validate(config.frequencyHz, config.background);
        rx.antenna.validate(config.frequencyHz, config.background);
        SimulationResult result{config, scene, tx, rx};
        const TraceOptions tracing{config.maxReflections, config.maxDistance};
        auto refinement = config.refinement;
        refinement.maxPathDistance = config.maxDistance;
        auto spreading = config.spreading;
        spreading.receiverTolerance = refinement.receiverTolerance;
        spreading.tMin = refinement.tMin;
        const auto &surfaces = scene.surfaces();

        // Sampled directions plus the exact direct direction; deduplication removes repeats.
        auto directions = launch_directions(config.rayCount);
        const Vec3 directOffset = rx.position - tx.position;
        if (directOffset.norm() > 0) directions.push_back(directOffset.normalized());
        // A planar reflection sequence has at most one specular solution; curved
        // sequences keep every discovery seed to find all stationary paths.
        std::set<std::vector<std::size_t>> convergedPlanar;
        for (const auto &direction: directions) {
            TracedRay ray(tx.position, direction);
            bool received = false;
            const auto traced = trace_ray(ray, surfaces, rx.position, config.receptionRadius, tracing,
                                          [&](const TracedRay &candidateRay) {
                                              received = true;
                                              std::vector<std::size_t> sequence;
                                              for (const auto &hit: candidateRay.reflections)
                                                  sequence.push_back(hit.surfaceIndex);
                                              const bool planar = std::none_of(
                                                  sequence.begin(), sequence.end(),
                                                  [&](auto i) { return curved(surfaces[i]); });
                                              if (planar && convergedPlanar.contains(sequence)) return;
                                              const Vec3 seed = sequence.empty() && directOffset.norm() > 0
                                                                    ? directOffset.normalized() : direction;
                                              auto refined = refine_path(tx.position, seed, rx.position, surfaces,
                                                                         sequence, refinement);
                                              if (planar && refined.status == RefinementStatus::Converged)
                                                  convergedPlanar.insert(sequence);
                                              result.candidates.push_back({std::move(refined), candidateRay.vertices});
                                          });
            const auto status = received ? TraceStatus::Received : traced.status;
            ++result.traceStatusCounts[static_cast<std::size_t>(status)];
        }

        std::vector<RefinementResult> trials;
        for (const auto &candidate: result.candidates) trials.push_back(candidate.refinement);
        FieldOptions fieldOptions{config.frequencyHz, config.background, config.sourceReference};
        for (auto index: deduplicate_paths(trials, config.deduplication)) {
            Path path;
            path.candidate = index;
            path.refinement = trials[index];
            const auto sequence = sequence_of(path.refinement.geometry);
            path.spreading = calculate_spreading(tx.position, path.refinement.launchDirection, rx.position, surfaces,
                                                 sequence, spreading);
            path.sourceField = tx.antenna.farfield(path.refinement.launchDirection, config.frequencyHz,
                                                   config.background) / spreading.referenceDistance;
            if (path.spreading.status == SpreadingStatus::Valid)
                path.field = reconstruct_field(path.refinement, path.spreading, path.sourceField, scene.materials(),
                                               fieldOptions);
            result.paths.push_back(std::move(path));
        }
        make_response(result);
        return result;
    }

    void set_receiver_antenna(SimulationResult &result, const Antenna &antenna) {
        result.receiver.antenna = antenna;
        make_response(result);
    }
} // namespace rt
