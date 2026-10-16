#pragma once
#include "material.hpp"
#include "refinement.hpp"
#include "spreading.hpp"

namespace rt {
    struct FieldOptions {
        double frequencyHz = 77e9;
        Material background; // homogeneous, lossless dielectric
        double sourceReference = 1.0; // common normalization for all launch directions
    };

    // Fresnel reflection at one hit; see Fresnel for the s/p conventions.
    struct Reflection {
        double incidenceAngle; // radians from the normal
        Complex s, p;
    };

    struct PathField {
        Vec3C field = Vec3C::Zero(); // receiver E (V/m)
        Vec3C normalizedField = Vec3C::Zero(); // field / sourceReference
        Vec3C transportedField = Vec3C::Zero(); // reflections only, before phase and spreading
        Vec3 arrivalDirection = Vec3::UnitX();
        std::vector<Reflection> reflections; // in path order
        double pathLength = 0.0, opticalPath = 0.0, delay = 0.0;
        double frequencyHz = 0.0, fieldFactor = 0.0;
        unsigned int causticCount = 0;
    };

    // Rebuild the receiver field along a converged path with valid spreading.
    // sourceField is transverse and given at spreading's reference distance.
    // Each reflection applies the Fresnel dyadic of surfaceMaterials[surface];
    // internal caustics contribute exp(-i pi count / 2).
    PathField reconstruct_field(const RefinementResult &path, const SpreadingResult &spreading,
                                const Vec3C &sourceField, const std::vector<Material> &surfaceMaterials,
                                const FieldOptions &options = {});

    struct ImpulseTap {
        double delay = 0.0; // seconds
        Complex coefficient = 0.0;
        Vec3C field = Vec3C::Zero(); // coherent sum of normalized receiver fields
        std::vector<std::size_t> paths; // contributing path indices
    };

    // Ideal complex-envelope response h(t) = sum coefficient_p delta(t - delay_p).
    struct ImpulseResponse {
        double frequencyHz = 0.0;
        std::vector<ImpulseTap> taps; // increasing delay
    };

    // Sort single-path taps by delay and coherently merge those within tolerance
    // of each group's earliest delay (no transitive chaining).
    ImpulseResponse group_taps(std::vector<ImpulseTap> paths, double frequencyHz, double delayTolerance = 1e-13);

    // Impulse response seen through a fixed unit (possibly complex) polarization
    // analyzer, using the Hermitian projection conj(analyzer) . field.
    ImpulseResponse calculate_impulse_response(const std::vector<PathField> &fields, const Vec3C &analyzer,
                                               double delayTolerance = 1e-13);

    // sum coefficient exp(+i 2 pi offset delay); at offset 0 the coherent carrier response.
    Complex frequency_response(const ImpulseResponse &response, double offsetHz = 0.0);
} // namespace rt
