#include "visualization.hpp"
#include <matplot/matplot.h>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

void plot_impulse(const ImpulseResponse &response, const std::filesystem::path &file) {
    std::vector<double> delays, magnitudes, phaseDelays, phases;
    double peak = 0;
    for (const auto &tap : response.taps_) peak = std::max(peak, std::abs(tap.coefficient_));
    for (const auto &tap : response.taps_) {
        delays.push_back(tap.delaySeconds_ * 1e9);
        magnitudes.push_back(std::abs(tap.coefficient_));
        if (std::abs(tap.coefficient_) > peak * 1e-12) {
            phaseDelays.push_back(tap.delaySeconds_ * 1e9);
            phases.push_back(std::arg(tap.coefficient_) * 180 / std::numbers::pi);
        }
    }
    auto figure = matplot::figure(true);
    figure->size(1000, 650);
    matplot::subplot(2, 1, 0);
    if (!delays.empty()) matplot::stem(delays, magnitudes);
    matplot::title("Impulse response at " + std::to_string(response.frequencyHz_ / 1e9) + " GHz");
    matplot::ylabel("Magnitude (common source reference)");
    matplot::grid(matplot::on);
    matplot::subplot(2, 1, 1);
    if (!phases.empty()) matplot::stem(phaseDelays, phases);
    matplot::xlabel("Absolute propagation delay (ns)");
    matplot::ylabel("Wrapped phase (degrees)");
    matplot::ylim({-180, 180});
    matplot::grid(matplot::on);
    if (!matplot::save(file.string())) throw std::runtime_error("Matplot++ could not save impulse plot");
}
