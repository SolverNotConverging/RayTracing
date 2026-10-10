#pragma once
#include "Antenna.hpp"
#include "Tracing.hpp"
#include "ImpulseResponse.hpp"

namespace rt {
    using ::frequency_response;
    using ::ImpulseResponse;
    using ::ImpulseTap;

    // Surfaces and the material assigned to each. Shapes and materials are
    // validated and copied when added; later edits to the caller's objects do
    // not affect the scene.
    class Scene {
    public:
        // Returns the surface index used by paths and reflection records.
        std::size_t add(Surface surface, Material material) {
            validate_surface(surface);
            material.validate();
            surfaces_.push_back(std::move(surface));
            materials_.push_back(std::move(material));
            return surfaces_.size() - 1;
        }

        const std::vector<Surface> &surfaces() const { return surfaces_; }

        // Material of each surface, indexed like surfaces().
        const std::vector<Material> &materials() const { return materials_; }

    private:
        std::vector<Surface> surfaces_;
        std::vector<Material> materials_;
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
        std::vector<std::size_t> traceStatusCounts = std::vector<std::size_t>(traceStatusCount, 0); // Indexed by TraceStatus
        std::vector<std::size_t> responseRayIndices; // Maps impulse pathIndices to rays.
    };

    SimulationResult solve(const Scene &scene, const Transmitter &tx,
                           const Receiver &rx, const SolverConfig &settings = {});

    // Reapply an Rx pattern/orientation to saved incident fields, without
    // retracing.
    void update_receiver(SimulationResult &result, const Antenna &antenna);
} // namespace rt
