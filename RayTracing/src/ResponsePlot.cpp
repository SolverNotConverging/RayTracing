#include "ResponsePlot.hpp"
#include "ResultIO.hpp"
#include "Visualization.hpp"
#include <algorithm>
#include <cmath>
#include <matplot/matplot.h>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace rt {
    void plot(const ImpulseResponse &response, const std::filesystem::path &png,
              bool showWindow) {
        if (!png.parent_path().empty())
            std::filesystem::create_directories(png.parent_path());
        double maximum = 0;
        for (const auto &tap: response.taps_)
            maximum = std::max(maximum, std::abs(tap.coefficient_));
        std::vector<double> delays, magnitudes, phaseDelays, phases;
        for (const auto &tap: response.taps_) {
            const double magnitude = std::abs(tap.coefficient_);
            delays.push_back(tap.delaySeconds_ * 1e9);
            magnitudes.push_back(magnitude);
            if (magnitude > maximum * 1e-12) {
                phaseDelays.push_back(delays.back());
                phases.push_back(std::arg(tap.coefficient_) * 180 / std::numbers::pi);
            }
        }
        using namespace matplot;
        auto f = figure(true);
        f->size(1200, 760);
        const double lastDelay =
                delays.empty() ? 1.0 : std::max(1.0, delays.back() * 1.05);
        auto magnitudeAxes = subplot(2, 1, 0);
        // Explicit margins reserve space for titles, ticks and the shared delay label.
        // Default subplot spacing can overlap text in the interactive Gnuplot window.
        magnitudeAxes->position({0.14f, 0.59f, 0.81f, 0.31f});
        magnitudeAxes->font_size(10);
        if (delays.empty())
            matplot::plot(std::vector<double>{0, lastDelay}, std::vector<double>{0, 0},
                          "b-");
        else
            stem(delays, magnitudes, "bo-")->line_width(1.5).marker_size(4);
        xlim({0, lastDelay});
        ylim({0, maximum > 0 ? maximum * 1.15 : 1.0});
        ylabel("Magnitude (common source reference)");
        std::ostringstream carrier;
        carrier << response.frequencyHz_ / 1e9;
        title(delays.empty()
                  ? "Impulse response: no valid field paths"
                  : "Impulse response magnitude (" + carrier.str() + " GHz)");
        grid(on);
        auto phaseAxes = subplot(2, 1, 1);
        phaseAxes->position({0.14f, 0.13f, 0.81f, 0.31f});
        phaseAxes->font_size(10);
        if (phaseDelays.empty())
            matplot::plot(std::vector<double>{0, lastDelay}, std::vector<double>{0, 0},
                          "w-");
        else
            stem(phaseDelays, phases, "ro-")->line_width(1.5).marker_size(4);
        xlim({0, lastDelay});
        ylim({-180, 180});
        xlabel("Absolute propagation delay (ns)");
        ylabel("Wrapped phase (degrees)");
        title(phaseDelays.empty()
                  ? "Phase undefined: no nonzero taps"
                  : "Impulse response phase");
        grid(on);
        if (!save(f, png.generic_string()))
            throw std::runtime_error("Cannot export impulse response plot");
        if (showWindow)
            show(f);
    }

    void plot_csv(const std::filesystem::path &csv,
                  const std::filesystem::path &png, bool showWindow) {
        plot(load_impulse_csv(csv), png, showWindow);
    }
} // namespace rt
void plot_impulse_response(const ImpulseResponse &response,
                           const std::vector<ReconstructedField> &,
                           const std::filesystem::path &directory,
                           bool showWindow) {
    rt::plot(response, directory / "impulse_response.png", showWindow);
}
