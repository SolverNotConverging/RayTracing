#include <raytracing/field.hpp>
#include <raytracing/core.hpp>
#include <raytracing/tracing.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace rt;

namespace {
constexpr double C = speedOfLight;
constexpr double PI = pi;

void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
void near(Complex actual, Complex expected, const char *message) {
    check(std::abs(actual - expected) < 1e-10, message);
}
template<class Action> void rejects(Action action) {
    bool rejected = false;
    try { action(); } catch (const std::invalid_argument &) { rejected = true; }
    check(rejected, "Invalid field/response input accepted");
}
PathField synthetic(double delay, Complex gain) {
    PathField field;
    field.delay = delay;
    field.frequencyHz = 77e9;
    field.normalizedField = gain * Vec3::UnitZ().cast<Complex>();
    return field;
}
}

int main() try {
    const Vec3 O = Vec3::Zero(), X = Vec3::UnitX(), Z = Vec3::UnitZ();
    const Vec3C z = Z.cast<Complex>();
    const auto path = refine_path(O, X, 3 * X, {}, {});
    const auto spreading = calculate_spreading(O, X, 3 * X, {}, {});
    FieldOptions options;
    options.frequencyHz = C / 12; // 3 m gives phase pi/2.
    const auto direct = reconstruct_field(path, spreading, z, {}, options);
    near(direct.field[2], Complex(0, 1.0 / 3), "Direct amplitude or phase incorrect");
    near(direct.delay, 3 / C, "Direct delay incorrect");
    const auto twice = reconstruct_field(path, spreading, 2 * z, {}, options);
    near(twice.field[2], 2.0 * direct.field[2], "Source field scaling incorrect");
    near(twice.normalizedField[2], 2.0 * direct.normalizedField[2], "Directional source amplitude was normalized away");
    options.background = Material::dielectric(2.25); // n = 1.5
    const auto medium = reconstruct_field(path, spreading, z, {}, options);
    near(medium.delay, 4.5 / C, "Refractive-index delay incorrect");
    near(medium.field[2], std::polar(1.0 / 3, 3 * PI / 4), "Optical path phase incorrect");
    near(medium.fieldFactor, direct.fieldFactor, "Homogeneous index altered geometric spreading");

    const Vec3 launch = Vec3(1, 1, 0).normalized();
    const std::vector<Surface> wall{Rectangle{X, Vec3::UnitY(), Z, 10, 10}};
    validate_surface(wall[0]);
    const auto reflectedPath = refine_path(O, launch, Vec3(0, 2, 0), wall, {0});
    const auto reflectedSpreading = calculate_spreading(O, launch, Vec3(0, 2, 0), wall, {0});
    options.background = {};
    const std::vector<Material> pec{Material::pec()};
    options.frequencyHz = C / (4 * std::sqrt(2.0)); // Total path gives pi phase.
    const auto reflected = reconstruct_field(reflectedPath, reflectedSpreading, z, pec, options);
    near(reflected.transportedField[2], -1, "PEC tangential sign incorrect");
    near(reflected.field[2], 1 / (2 * std::sqrt(2.0)), "Reflection and propagation phase were not combined correctly");
    near(reflected.delay, 2 * std::sqrt(2.0) / C, "Reflected delay uses sphere entry instead of full path");
    const Vec3 u = launch.unitOrthogonal(), v = launch.cross(u).normalized();
    const Vec3C circular = (u.cast<Complex>() + Complex(0, 1) * v.cast<Complex>()) / std::sqrt(2.0);
    const auto polarized = reconstruct_field(reflectedPath, reflectedSpreading, circular, pec, options);
    const Vec3C normal = X.cast<Complex>();
    const Vec3C expected = -circular + 2.0 * normal.dot(circular) * normal;
    check((polarized.transportedField - expected).norm() < 1e-10, "Complex PEC polarization transport incorrect");
    check(std::abs(polarized.arrivalDirection.cast<Complex>().dot(polarized.field)) < 1e-10,
        "Reconstructed field is not transverse");

    const auto single = calculate_impulse_response({direct}, z);
    near(single.taps[0].coefficient, direct.field[2], "Carrier phase missing from impulse coefficient");
    near(frequency_response(single), direct.field[2], "Carrier response differs from coherent tap sum");
    near(frequency_response(single, C / 12), -1.0 / 3, "Envelope frequency convention inconsistent with propagation");
    const auto interference = calculate_impulse_response({synthetic(2e-9, 0.5), synthetic(1e-9, 1),
        synthetic(1e-9, Complex(0, 1))}, z);
    check(interference.taps.size() == 2 && interference.taps[0].paths.size() == 2,
        "Equal-delay fields were not coherently grouped and sorted");
    near(interference.taps[0].coefficient, Complex(1, 1), "Equal-delay complex addition incorrect");
    near(frequency_response(interference), Complex(1.5, 1), "Multipath carrier response incorrect");
    const auto cancelled = calculate_impulse_response({synthetic(1e-9, 1), synthetic(1e-9, -1)}, z);
    near(cancelled.taps[0].coefficient, 0, "Destructive interference lost");
    const auto groups = calculate_impulse_response({synthetic(1e-9, 1), synthetic(1e-9 + 0.75e-13, 1),
        synthetic(1e-9 + 1.5e-13, 1)}, z);
    check(groups.taps.size() == 2, "Delay grouping chained transitively");
    const Vec3C analyzer = (Vec3::UnitY().cast<Complex>() + Complex(0, 1) * z) / std::sqrt(2.0);
    auto circularTap = synthetic(1e-9, 0);
    circularTap.normalizedField = Complex(0, 1) * analyzer;
    near(calculate_impulse_response({circularTap}, analyzer).taps[0].coefficient, Complex(0, 1),
        "Receiver analyzer did not use Hermitian projection");
    const auto empty = calculate_impulse_response({}, z);
    check(empty.taps.empty(), "Empty field list generated taps");
    near(frequency_response(empty), 0, "Empty channel is not zero");

    auto invalidPath = path;
    invalidPath.status = RefinementStatus::UnresolvedCorner;
    rejects([&] { reconstruct_field(invalidPath, spreading, z, {}, options); });
    auto caustic = spreading;
    caustic.status = SpreadingStatus::Caustic;
    rejects([&] { reconstruct_field(path, caustic, z, {}, options); });
    std::cout << "All field reconstruction and impulse response checks passed\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
