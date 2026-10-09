#include "ImpulseResponse.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <numbers>
#include <stdexcept>

ImpulseResponse calculate_impulse_response(const std::vector<ReconstructedField> &fields,
                                           const Vec3C &receiverPolarization, double delayToleranceSeconds) {
    if (!receiverPolarization.allFinite() || std::abs(receiverPolarization.norm() - 1.0) > 1e-10 ||
        !std::isfinite(delayToleranceSeconds) || delayToleranceSeconds < 0)
        throw std::invalid_argument(
            "Impulse response requires unit receiver polarization and nonnegative delay tolerance");
    ImpulseResponse response;
    response.receiverPolarization_ = receiverPolarization;
    response.delayToleranceSeconds_ = delayToleranceSeconds;
    if (!fields.empty()) response.frequencyHz_ = fields.front().frequencyHz_;
    std::vector<std::size_t> indices(fields.size());
    std::iota(indices.begin(), indices.end(), 0);
    for (const auto &field: fields) {
        if (!std::isfinite(field.delaySeconds_) || field.delaySeconds_ < 0 ||
            !field.normalizedReceiverField_.allFinite() || !std::isfinite(field.frequencyHz_) || field.frequencyHz_ <= 0
            ||
            std::abs(field.frequencyHz_ - response.frequencyHz_) > 1e-12 * response.frequencyHz_)
            throw std::invalid_argument("Impulse response fields must be finite and use the same carrier frequency");
    }
    std::stable_sort(indices.begin(), indices.end(), [&](std::size_t a, std::size_t b) {
        return fields[a].delaySeconds_ < fields[b].delaySeconds_;
    });
    for (std::size_t index: indices) {
        const auto &field = fields[index];
        if (response.taps_.empty() || field.delaySeconds_ - response.taps_.back().delaySeconds_ > delayToleranceSeconds)
            response.taps_.push_back({field.delaySeconds_, {}, Vec3C::Zero(), {}});
        auto &tap = response.taps_.back();
        tap.normalizedField_ += field.normalizedReceiverField_;
        tap.pathIndices_.push_back(index);
    }
    for (auto &tap: response.taps_) {
        tap.coefficient_ = receiverPolarization.dot(tap.normalizedField_);
        if (!tap.normalizedField_.allFinite() || !std::isfinite(tap.coefficient_.real()) || !std::isfinite(
                tap.coefficient_.imag()))
            throw std::invalid_argument("Coherent impulse tap overflow");
    }
    return response;
}

Complex frequency_response(const ImpulseResponse &response, double offsetHz) {
    if (!std::isfinite(offsetHz)) throw std::invalid_argument("Frequency offset must be finite");
    Complex sum = 0.0;
    for (const auto &tap: response.taps_) {
        const double cycles = offsetHz * tap.delaySeconds_;
        if (!std::isfinite(cycles)) throw std::invalid_argument("Frequency response phase overflow");
        sum += tap.coefficient_ * std::polar(1.0, 2.0 * std::numbers::pi * std::remainder(cycles, 1.0));
    }
    return sum;
}
