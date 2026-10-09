#pragma once
#include "Antenna.hpp"
#include "EMRay.hpp"
#include "ImpulseResponse.hpp"

namespace rt {
    using ::frequency_response;
    using ::ImpulseResponse;
    using ::ImpulseTap;

    class Scene {
    public:
        void add(Surface surface) {
            validate_surface(surface);
            surfaces_.push_back(std::move(surface));
        }

        const std::vector<Surface> &surfaces() const { return surfaces_; }

    private:
        std::vector<Surface> surfaces_;
    };

    struct Transmitter {
        Vec3 position = Vec3::Zero();
        Antenna antenna;
    };

    struct Receiver {
        Vec3 position = Vec3::Zero();
        Antenna antenna;
    };

    struct SolverConfig {
        double frequencyHz = 77e9;
        std::size_t rayCount = 2000;
        int maxReflections = 8;
        double maxDistance = 20.0, receptionRadius = 0.12;
        double commonSourceReference = 1.0;
        double delayToleranceSeconds = 1e-13;
        Medium medium;
        RefinementOptions refinement;
        SpreadingOptions spreading;
        DeduplicationOptions deduplication;
    };

    struct Candidate {
        RefinementResult refinement;
        std::vector<Vec3> coarseVertices;
    };

    struct SolvedRay {
        std::size_t candidateIndex = 0;
        RefinementResult refinement;
        SpreadingResult spreading;
        std::optional<ReconstructedField> field;
        Vec3C sourceField = Vec3C::Zero();
        Complex receivedCoefficient = 0;
    };

    struct SimulationResult {
        SolverConfig settings;
        Scene scene;
        Transmitter transmitter;
        Receiver receiver;
        std::vector<SolvedRay> rays;
        std::vector<Candidate> candidates;
        ImpulseResponse impulseResponse;
        std::size_t launchedRays = 0;
        std::vector<std::size_t> traceStatusCounts = std::vector<std::size_t>(6, 0);
        std::vector<std::size_t>
        responseRayIndices; // Maps impulse pathIndices to rays.
    };

    SimulationResult solve(const Scene &scene, const Transmitter &tx,
                           const Receiver &rx, const SolverConfig &settings = {});

    // Reapply an Rx pattern/orientation to saved incident fields, without
    // retracing.
    void update_receiver(SimulationResult &result, const Antenna &antenna);
} // namespace rt
