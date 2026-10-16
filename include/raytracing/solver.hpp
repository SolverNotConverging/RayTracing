#pragma once
#include "antenna.hpp"
#include "field.hpp"
#include "tracing.hpp"

// Point-source multipath solver: discover reflection sequences with coarse
// rays, refine exact specular paths, then rebuild coherent receiver fields.
namespace rt {
    // Reflecting surfaces, each the boundary of a half-space of its material.
    class Scene {
    public:
        // Validate and copy; returns the surface index used by paths.
        std::size_t add(Surface surface, Material material);
        const std::vector<Surface> &surfaces() const { return surfaces_; }
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
        double maxDistance = 20.0;
        double receptionRadius = 0.12;
        double sourceReference = 1.0;
        double delayTolerance = 1e-13; // seconds
        Material background;
        RefinementOptions refinement;
        SpreadingOptions spreading;
        DeduplicationOptions deduplication;
    };

    struct Candidate {
        RefinementResult refinement;
        std::vector<Vec3> coarseVertices;
    };

    struct Path {
        std::size_t candidate = 0;
        RefinementResult refinement;
        SpreadingResult spreading;
        std::optional<PathField> field; // absent when spreading is singular
        Vec3C sourceField = Vec3C::Zero();
        Complex coefficient = 0.0; // receive-antenna response

        std::vector<Vec3> vertices() const; // Tx, reflections, Rx
    };

    struct SimulationResult {
        SolverConfig config;
        Scene scene;
        Transmitter transmitter;
        Receiver receiver;
        std::vector<Path> paths;
        // Discovery trials; converged planar sequences keep one representative.
        std::vector<Candidate> candidates;
        ImpulseResponse impulseResponse; // tap.paths index into paths
        std::vector<std::size_t> traceStatusCounts = std::vector<std::size_t>(traceStatusCount, 0);
    };

    SimulationResult solve(const Scene &scene, const Transmitter &tx, const Receiver &rx,
                           const SolverConfig &config = {});

    // Reapply a receive antenna to the stored incident fields without retracing.
    void set_receiver_antenna(SimulationResult &result, const Antenna &antenna);
} // namespace rt
