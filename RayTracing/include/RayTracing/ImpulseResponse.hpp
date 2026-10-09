#pragma once

#include "FieldReconstruction.hpp"
#include <string>

struct ImpulseTap {
    double delaySeconds_ = 0.0;
    Complex coefficient_ = 0.0; // Receiver polarization projection, normalized source amplitude.
    Vec3C normalizedField_ = Vec3C::Zero();
    std::vector<std::size_t> pathIndices_; // Indices into the supplied field vector.
};

struct ImpulseResponse {
    double frequencyHz_ = 0.0;
    Vec3C receiverPolarization_ = Vec3C::Zero(); // Unit Hermitian analyzer in global coordinates.
    std::vector<ImpulseTap> taps_; // Increasing absolute propagation delay.
    double delayToleranceSeconds_ = 1e-13;
    std::string receiverModel_ = "fixed Hermitian polarization analyzer";
};

// Ideal complex envelope response h(t) = sum coefficient_p delta(t-delay_p).
// Coherently sums delays within tolerance of each group's earliest delay (no
// transitive grouping). Tolerance is numerical, not a receiver bandwidth/bin size.
// Fixed receiver polarization must be unit length, and may be complex.
ImpulseResponse calculate_impulse_response(const std::vector<ReconstructedField> &fields,
                                           const Vec3C &receiverPolarization, double delayToleranceSeconds = 1e-13);

// With the exp(-i omega t) convention, the envelope frequency response is
// sum coefficient_p exp(+i 2pi offsetHz delay_p). At offset=0 it is the coherent
// normalized receiver field at the carrier, not a sum of magnitudes.
Complex frequency_response(const ImpulseResponse &response, double offsetHz = 0.0);
