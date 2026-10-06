#include <Eigen/Dense>
#include <matplot/matplot.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

// ============================================================
// Basic types
// ============================================================

using Vec3 = Eigen::Vector3d;

constexpr double PI = 3.14159265358979323846;

// ============================================================
// Simulation data
// ============================================================

struct SimulationConfig
{
    double sample_rate = 48000.0;      // Hz
    double duration = 1.0;             // seconds
    double sound_speed = 343.0;        // m/s
};

struct Source
{
    Vec3 position;
};

struct Receiver
{
    Vec3 position;

    // We won't use this until ray tracing starts,
    // but define it now because the receiver will
    // eventually be represented by a small sphere.
    double radius = 0.10;              // metres
};

struct Room
{
    Vec3 min_corner;
    Vec3 max_corner;
};

// ============================================================
// Ray representation
// ============================================================

struct Ray
{
    Vec3 origin;
    Vec3 direction;

    // Distance travelled by the ray before its current origin.
    double travelled_distance = 0.0;

    // Pressure-like multiplicative weight.
    double amplitude = 1.0;
};

// ============================================================
// IR helper
// ============================================================

void depositImpulse(
    std::vector<double>& ir,
    double arrival_time,
    double amplitude,
    double sample_rate)
{
    // Continuous-valued sample position.
    const double sample_position = arrival_time * sample_rate;

    if (sample_position < 0.0)
        return;

    const std::size_t i =
        static_cast<std::size_t>(std::floor(sample_position));

    const double fraction =
        sample_position - static_cast<double>(i);

    if (i >= ir.size())
        return;

    // Linear fractional-delay interpolation.
    ir[i] += amplitude * (1.0 - fraction);

    if (i + 1 < ir.size())
        ir[i + 1] += amplitude * fraction;
}

// ============================================================
// Main
// ============================================================

int main()
{
    SimulationConfig config;

    Room room{
        .min_corner = Vec3{0.0, 0.0, 0.0},
        .max_corner = Vec3{6.0, 5.0, 3.0}
    };

    Source source{
        .position = Vec3{1.0, 1.0, 1.5}
    };

    Receiver receiver{
        .position = Vec3{4.0, 3.0, 1.5},
        .radius = 0.10
    };

    // --------------------------------------------------------
    // Allocate the impulse response
    // --------------------------------------------------------

    const std::size_t num_samples =
        static_cast<std::size_t>(
            config.duration * config.sample_rate);

    std::vector<double> ir(num_samples, 0.0);

    // --------------------------------------------------------
    // Direct path
    // --------------------------------------------------------

    const Vec3 source_to_receiver =
        receiver.position - source.position;

    const double direct_distance =
        source_to_receiver.norm();

    const double arrival_time =
        direct_distance / config.sound_speed;

    const double amplitude =
        1.0 / (4.0 * PI * direct_distance);

    depositImpulse(
        ir,
        arrival_time,
        amplitude,
        config.sample_rate
    );

    // --------------------------------------------------------
    // Diagnostics
    // --------------------------------------------------------

    std::cout << "Direct path distance : "
              << direct_distance << " m\n";

    std::cout << "Direct arrival time  : "
              << arrival_time * 1000.0 << " ms\n";

    std::cout << "Continuous sample    : "
              << arrival_time * config.sample_rate << '\n';

    std::cout << "Direct amplitude     : "
              << amplitude << '\n';

    // --------------------------------------------------------
    // Time axis for plotting
    // --------------------------------------------------------

    std::vector<double> time_ms(num_samples);

    for (std::size_t i = 0; i < num_samples; ++i)
    {
        time_ms[i] =
            1000.0 *
            static_cast<double>(i) /
            config.sample_rate;
    }

    // --------------------------------------------------------
    // Plot
    // --------------------------------------------------------

    matplot::plot(time_ms, ir);

    matplot::xlabel("Time [ms]");
    matplot::ylabel("Amplitude");
    matplot::title("Room impulse response - direct path only");

    matplot::show();

    return 0;
}