#include "FieldReconstruction.hpp"
#include "ImpulseResponse.hpp"
#include "Constants.hpp"
#include "Tracing.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
constexpr double C = rt::constants::speedOfLight;
constexpr double PI = rt::constants::pi;

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
ReconstructedField synthetic(double delay, Complex gain) {
    ReconstructedField field;
    field.delaySeconds_ = delay;
    field.frequencyHz_ = 77e9;
    field.normalizedReceiverField_ = gain * Vec3::UnitZ().cast<Complex>();
    return field;
}
}

int main() try {
    const Vec3 O = Vec3::Zero(), X = Vec3::UnitX(), Z = Vec3::UnitZ();
    const Vec3C z = Z.cast<Complex>();
    const auto path = refine_path(O, X, 3 * X, {}, {});
    const auto spreading = calculate_spreading(O, X, 3 * X, {}, {});
    FieldReconstructionOptions options;
    options.frequencyHz_ = C / 12; // 3 m gives phase pi/2.
    const auto direct = reconstruct_field(path, spreading, z, {}, options);
    near(direct.receiverField_[2], Complex(0, 1.0 / 3), "Direct amplitude or phase incorrect");
    near(direct.delaySeconds_, 3 / C, "Direct delay incorrect");
    const auto twice = reconstruct_field(path, spreading, 2 * z, {}, options);
    near(twice.receiverField_[2], 2.0 * direct.receiverField_[2], "Source field scaling incorrect");
    near(twice.normalizedReceiverField_[2], 2.0 * direct.normalizedReceiverField_[2], "Directional source amplitude was normalized away");
    options.medium_.relativePermittivity = 2.25; // n = 1.5
    const auto medium = reconstruct_field(path, spreading, z, {}, options);
    near(medium.delaySeconds_, 4.5 / C, "Refractive-index delay incorrect");
    near(medium.receiverField_[2], std::polar(1.0 / 3, 3 * PI / 4), "Optical path phase incorrect");
    near(medium.fieldFactor_, direct.fieldFactor_, "Homogeneous index altered geometric spreading");

    const Vec3 launch = Vec3(1, 1, 0).normalized();
    const std::vector<Surface> wall{Rectangle{X, Vec3::UnitY(), Z, 10, 10}};
    validate_surface(wall[0]);
    const auto reflectedPath = refine_path(O, launch, Vec3(0, 2, 0), wall, {0});
    const auto reflectedSpreading = calculate_spreading(O, launch, Vec3(0, 2, 0), wall, {0});
    options.medium_ = {};
    const std::vector<rt::Material> pec{rt::Material::pec()};
    options.frequencyHz_ = C / (4 * std::sqrt(2.0)); // Total path gives pi phase.
    const auto reflected = reconstruct_field(reflectedPath, reflectedSpreading, z, pec, options);
    near(reflected.transportedReferenceField_[2], -1, "PEC tangential sign incorrect");
    near(reflected.receiverField_[2], 1 / (2 * std::sqrt(2.0)), "Reflection and propagation phase were not combined correctly");
    near(reflected.delaySeconds_, 2 * std::sqrt(2.0) / C, "Reflected delay uses sphere entry instead of full path");
    const Vec3 u = launch.unitOrthogonal(), v = launch.cross(u).normalized();
    const Vec3C circular = (u.cast<Complex>() + Complex(0, 1) * v.cast<Complex>()) / std::sqrt(2.0);
    const auto polarized = reconstruct_field(reflectedPath, reflectedSpreading, circular, pec, options);
    const Vec3C normal = X.cast<Complex>();
    const Vec3C expected = -circular + 2.0 * normal.dot(circular) * normal;
    check((polarized.transportedReferenceField_ - expected).norm() < 1e-10, "Complex PEC polarization transport incorrect");
    check(std::abs(polarized.arrivalDirection_.cast<Complex>().dot(polarized.receiverField_)) < 1e-10,
        "Reconstructed field is not transverse");

    const auto single = calculate_impulse_response({direct}, z);
    near(single.taps_[0].coefficient_, direct.receiverField_[2], "Carrier phase missing from impulse coefficient");
    near(frequency_response(single), direct.receiverField_[2], "Carrier response differs from coherent tap sum");
    near(frequency_response(single, C / 12), -1.0 / 3, "Envelope frequency convention inconsistent with propagation");
    const auto interference = calculate_impulse_response({synthetic(2e-9, 0.5), synthetic(1e-9, 1),
        synthetic(1e-9, Complex(0, 1))}, z);
    check(interference.taps_.size() == 2 && interference.taps_[0].pathIndices_.size() == 2,
        "Equal-delay fields were not coherently grouped and sorted");
    near(interference.taps_[0].coefficient_, Complex(1, 1), "Equal-delay complex addition incorrect");
    near(frequency_response(interference), Complex(1.5, 1), "Multipath carrier response incorrect");
    const auto cancelled = calculate_impulse_response({synthetic(1e-9, 1), synthetic(1e-9, -1)}, z);
    near(cancelled.taps_[0].coefficient_, 0, "Destructive interference lost");
    const auto groups = calculate_impulse_response({synthetic(1e-9, 1), synthetic(1e-9 + 0.75e-13, 1),
        synthetic(1e-9 + 1.5e-13, 1)}, z);
    check(groups.taps_.size() == 2, "Delay grouping chained transitively");
    const Vec3C analyzer = (Vec3::UnitY().cast<Complex>() + Complex(0, 1) * z) / std::sqrt(2.0);
    auto circularTap = synthetic(1e-9, 0);
    circularTap.normalizedReceiverField_ = Complex(0, 1) * analyzer;
    near(calculate_impulse_response({circularTap}, analyzer).taps_[0].coefficient_, Complex(0, 1),
        "Receiver analyzer did not use Hermitian projection");
    const auto empty = calculate_impulse_response({}, z);
    check(empty.taps_.empty(), "Empty field list generated taps");
    near(frequency_response(empty), 0, "Empty channel is not zero");

    rejects([&] { reconstruct_field(path, spreading, X.cast<Complex>(), {}, options); });
    auto invalidPath = path;
    invalidPath.status_ = RefinementStatus::UnresolvedCorner;
    rejects([&] { reconstruct_field(invalidPath, spreading, z, {}, options); });
    auto caustic = spreading;
    caustic.status_ = SpreadingStatus::Caustic;
    rejects([&] { reconstruct_field(path, caustic, z, {}, options); });
    rejects([&] { reconstruct_field(reflectedPath, reflectedSpreading, z, {}, options); }); // No material
    rejects([&] { calculate_impulse_response({direct}, 2 * z); });
    auto wrongFrequency = direct;
    wrongFrequency.frequencyHz_ *= 2;
    rejects([&] { calculate_impulse_response({direct, wrongFrequency}, z); });
    rejects([&] { calculate_impulse_response({}, z, -1); });
    std::cout << "All field reconstruction and impulse response checks passed\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
