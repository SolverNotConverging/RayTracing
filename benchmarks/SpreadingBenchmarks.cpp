#include "Spreading.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

namespace {
    struct Case {
        Vec3 tx, direction, rx;
        std::vector<Surface> surfaces;
        std::vector<std::size_t> sequence;
        Eigen::Vector2d expectedLengths;
    };

    void check(bool condition, const std::string &message) {
        if (!condition) throw std::runtime_error(message);
    }

    double relative(double actual, double expected) {
        return std::abs(actual - expected) / std::abs(expected);
    }

    Case curved(double a, double b, double radius, double angle, bool cylinder = false, bool concave = false) {
        const double c = std::cos(angle), s = std::sin(angle);
        const Vec3 incoming{concave ? c : -c, s, 0};
        const Vec3 outgoing{-incoming.x(), s, 0};
        const double sign = concave ? -1.0 : 1.0;
        Case scene{
            -a * incoming, incoming, b * outgoing, {}, {0},
            {
                a + b + sign * 2 * a * b / (radius * c),
                cylinder ? a + b : a + b + sign * 2 * a * b * c / radius
            }
        };
        if (cylinder) scene.surfaces.push_back(Cylinder{Vec3(-radius, 0, 0), Vec3::UnitZ(), radius, 20, true});
        else scene.surfaces.push_back(Sphere{Vec3(-radius, 0, 0), radius});
        return scene;
    }

    Case direct(double distance) {
        return {Vec3::Zero(), Vec3::UnitX(), Vec3(distance, 0, 0), {}, {}, {distance, distance}};
    }

    Case plane(double b) {
        return {
            Vec3(2, 0, 0), -Vec3::UnitX(), Vec3(b, 0, 0),
            {Rectangle{Vec3::Zero(), Vec3::UnitY(), Vec3::UnitZ(), 20, 20}}, {0}, {2 + b, 2 + b}
        };
    }

    Case two_planes(double b) {
        const Vec3 d = Vec3(1, 0.25, 0).normalized();
        const Vec3 secondHit(-1, 0.75, 0);
        const double length = 3 / d.x() + b;
        return {
            Vec3::Zero(), d, secondHit + b * d,
            {
                Rectangle{Vec3(1, 0, 0), Vec3::UnitY(), Vec3::UnitZ(), 20, 20},
                Rectangle{Vec3(-1, 0, 0), Vec3::UnitY(), Vec3::UnitZ(), 20, 20}
            },
            {0, 1}, {length, length}
        };
    }

    Case two_spheres(double b) {
        // Unfolded paraxial transfer: P(b) M(R2) P(d) M(R1) P(a).
        // Propagation changes height by length*slope; a convex reflection
        // changes slope by 2*height/R. Start with height=0, slope=1.
        constexpr double a = 1, d = 3, r1 = 1, r2 = 0.7;
        const double firstSlope = 1 + 2 * a / r1;
        const double secondHeight = a + d * firstSlope;
        const double secondSlope = firstSlope + 2 * secondHeight / r2;
        const double length = secondHeight + b * secondSlope;
        return {
            Vec3(-a, 0, 0), Vec3::UnitX(), Vec3(-d + b, 0, 0),
            {Sphere{Vec3(r1, 0, 0), r1}, Sphere{Vec3(-d - r2, 0, 0), r2}},
            {0, 1}, {length, length}
        };
    }

    Case transformed(Case scene) {
        const Eigen::Matrix3d rotation = Eigen::AngleAxisd(0.73, Vec3(1, 2, 3).normalized()).toRotationMatrix();
        const Vec3 shift(3, -2, 1);
        scene.tx = rotation * scene.tx + shift;
        scene.rx = rotation * scene.rx + shift;
        scene.direction = rotation * scene.direction;
        for (auto &surface: scene.surfaces) {
            std::visit([&](auto &shape) {
                using T = std::decay_t<decltype(shape)>;
                if constexpr (std::is_same_v<T, Triangle>) {
                    shape.a_ = rotation * shape.a_ + shift;
                    shape.b_ = rotation * shape.b_ + shift;
                    shape.c_ = rotation * shape.c_ + shift;
                } else {
                    shape.center_ = rotation * shape.center_ + shift;
                    if constexpr (std::is_same_v<T, Rectangle>) {
                        shape.u_ = rotation * shape.u_;
                        shape.v_ = rotation * shape.v_;
                    } else if constexpr (std::is_same_v<T, Disk>) shape.normal_ = rotation * shape.normal_;
                    else if constexpr (std::is_same_v<T, Cylinder>) shape.axis_ = rotation * shape.axis_;
                }
            }, surface);
        }
        return scene;
    }

    SpreadingResult evaluate(const Case &scene, const SpreadingOptions &options = {}) {
        for (const auto &surface: scene.surfaces) validate_surface(surface);
        return calculate_spreading(scene.tx, scene.direction, scene.rx, scene.surfaces, scene.sequence, options);
    }

    double validate(const Case &scene, const SpreadingResult &result, double tolerance = 2e-6) {
        check(result.status_ == SpreadingStatus::Valid && result.fieldFactor_,
              std::string("Unexpected spreading status: ") + spreading_status_name(result.status_));
        Eigen::Vector2d lengths = scene.expectedLengths.cwiseAbs();
        if (lengths[0] < lengths[1]) std::swap(lengths[0], lengths[1]);
        double error = std::max(relative(result.singularValues_[0], lengths[0]),
                                relative(result.singularValues_[1], lengths[1]));
        error = std::max(error, relative(result.areaPerSolidAngle_, lengths.prod()));
        error = std::max(error, relative(*result.fieldFactor_, 1.0 / std::sqrt(lengths.prod())));
        check(error < tolerance, "Analytical spreading mismatch: " + std::to_string(error));
        return error;
    }
}

int main(int argc, char **argv) try {
    if (argc > 2) throw std::invalid_argument("Usage: SpreadingBenchmarks [output_directory]");
    const std::filesystem::path output = argc == 2 ? argv[1] : "benchmarks/results";
    std::filesystem::create_directories(output);
    std::ofstream csv(output / "spreading.csv");
    check(bool(csv), "Cannot open benchmark CSV");
    csv << std::setprecision(17)
            <<
            "case,distance_m,analytical_sigma_max_m,analytical_sigma_min_m,numerical_sigma_max_m,numerical_sigma_min_m,analytical_area_m2,numerical_area_m2,analytical_field_factor,numerical_field_factor,max_relative_error\n";
    int checks = 0;
    double maximumError = 0;
    const std::vector<std::string> names{
        "direct", "plane", "sphere_normal", "cylinder_normal",
        "sphere_oblique_45deg", "sphere_oblique_70deg", "cylinder_oblique_45deg", "sphere_concave", "two_planes",
        "two_spheres"
    };
    for (const auto &name: names) {
        double caseError = 0;
        for (int i = 0; i <= 60; ++i) {
            const double b = name == "two_planes"
                                 ? 0.1 + i * 0.02
                                 : (name == "sphere_concave" || name == "two_spheres")
                                       ? 0.1 + i * 0.04
                                       : 0.2 + i * 0.13;
            Case scene = direct(b);
            if (name == "plane") scene = plane(b);
            else if (name == "sphere_normal") scene = curved(2, b, 1, 0);
            else if (name == "cylinder_normal") scene = curved(2, b, 1, 0, true);
            else if (name == "sphere_oblique_45deg") scene = curved(2, b, 1, std::numbers::pi / 4);
            else if (name == "sphere_oblique_70deg") scene = curved(2, b, 1, 70 * std::numbers::pi / 180);
            else if (name == "cylinder_oblique_45deg") scene = curved(2, b, 1, std::numbers::pi / 4, true);
            else if (name == "sphere_concave") scene = curved(3, b, 2, 0, false, true);
            else if (name == "two_planes") scene = two_planes(b);
            else if (name == "two_spheres") scene = two_spheres(b);
            if (scene.expectedLengths.cwiseAbs().minCoeff() < 1e-6) {
                csv << name << ',' << b << ",nan,nan,nan,nan,nan,nan,nan,nan,nan\n";
                continue; // Separate caustic test below.
            }
            const auto result = evaluate(scene);
            const double error = validate(scene, result);
            const auto moved = evaluate(transformed(scene));
            const double movedError = validate(scene, moved);
            check(relative(result.areaPerSolidAngle_, moved.areaPerSolidAngle_) < 2e-6,
                  "Rigid transformation changed spreading");
            checks += 2;
            caseError = std::max({caseError, error, movedError});
            Eigen::Vector2d lengths = scene.expectedLengths.cwiseAbs();
            if (lengths[0] < lengths[1]) std::swap(lengths[0], lengths[1]);
            const double expectedFactor = 1 / std::sqrt(lengths.prod());
            csv << name << ',' << b << ',' << lengths[0] << ',' << lengths[1] << ','
                    << result.singularValues_[0] << ',' << result.singularValues_[1] << ','
                    << lengths.prod() << ',' << result.areaPerSolidAngle_ << ',' << expectedFactor << ','
                    << *result.fieldFactor_ << ',' << error << '\n';
        }
        maximumError = std::max(maximumError, caseError);
        std::cout << name << ": maximum relative error " << caseError << '\n';
    }

    // Near and at a concave sphere's axial focus: b = a R / (2a - R).
    const auto focus = evaluate(curved(3, 1.5, 2, 0, false, true));
    check(focus.status_ == SpreadingStatus::Caustic && !focus.fieldFactor_, "Caustic must not supply an amplitude");
    const auto cylinderFocus = evaluate(curved(3, 1.5, 2, 0, true, true));
    check(cylinderFocus.status_ == SpreadingStatus::Caustic && !cylinderFocus.fieldFactor_,
          "One-axis caustic must not supply an amplitude");
    SpreadingOptions reference;
    reference.referenceDistance_ = 2;
    const auto scaled = evaluate(direct(5), reference);
    check(scaled.fieldFactor_ && std::abs(*scaled.fieldFactor_ - 0.4) < 1e-10, "Reference-distance scaling failed");
    auto miss = direct(5);
    miss.rx.y() = 0.01;
    check(evaluate(miss).status_ == SpreadingStatus::InvalidCentralPath, "Unrefined ray accepted");
    auto blocked = direct(5);
    blocked.surfaces.push_back(Sphere{Vec3(2, 0, 0), 0.1});
    check(evaluate(blocked).status_ == SpreadingStatus::InvalidCentralPath, "Blocked central ray accepted");
    // A central ray on an aperture edge has no two-sided smooth neighbourhood.
    Case edge{
        Vec3(1, 0, 0), -Vec3::UnitX(), Vec3(1, 0, 0),
        {Rectangle{Vec3(0, 1, 0), Vec3::UnitY(), Vec3::UnitZ(), 1, 1}}, {0}, {2, 2}
    };
    check(evaluate(edge).status_ == SpreadingStatus::InvalidPerturbation, "Aperture-edge derivative accepted");
    bool rejected = false;
    try {
        reference.angularStep_ = 0;
        evaluate(direct(5), reference);
    } catch (const std::invalid_argument &) { rejected = true; }
    check(rejected, "Invalid options accepted");
    checks += 7;

    SpreadingOptions strict;
    strict.maxStepHalvings_ = 0;
    strict.relativeDerivativeTolerance_ = 1e-14;
    strict.absoluteDerivativeTolerance_ = 1e-14;
    const auto unstable = evaluate(curved(2, 3, 1, std::numbers::pi / 3), strict);
    check(unstable.status_ == SpreadingStatus::UnstableDerivative && !unstable.fieldFactor_,
          "Unstable derivative supplied an amplitude");
    auto neighbouringBlocker = direct(5);
    neighbouringBlocker.surfaces.push_back(Sphere{Vec3(2, 0.0002, 0), 0.00005});
    const auto adapted = evaluate(neighbouringBlocker);
    validate(neighbouringBlocker, adapted);
    check(adapted.angularStep_ < SpreadingOptions{}.angularStep_ / 2,
          "Invalid neighbouring ray did not trigger step reduction");
    checks += 2;

    // Measure truncation error without adaptive halving masking the h^2 trend.
    std::ofstream convergence(output / "convergence.csv");
    check(bool(convergence), "Cannot open convergence CSV");
    convergence << std::setprecision(17) << "actual_angular_step,relative_area_error\n";
    std::vector<double> convergenceErrors;
    const auto scene = curved(2, 3, 1, std::numbers::pi / 3);
    for (double h: {0.02, 0.01, 0.005, 0.002, 0.001, 0.0002, 0.0001, 0.00002, 0.00001, 0.000002}) {
        SpreadingOptions options;
        options.angularStep_ = h;
        options.maxStepHalvings_ = 0;
        options.relativeDerivativeTolerance_ = 1.0;
        const auto result = evaluate(scene, options);
        check(result.status_ == SpreadingStatus::Valid, "Convergence sample failed");
        const double error = relative(result.areaPerSolidAngle_, scene.expectedLengths.prod());
        convergenceErrors.push_back(std::max(error, 1e-16));
        convergence << result.angularStep_ << ',' << error << '\n';
    }
    const double errorRatio = convergenceErrors[0] / convergenceErrors[1];
    check(errorRatio > 3.5 && errorRatio < 4.5, "Central differences did not show second-order convergence");
    ++checks;
    check(bool(csv) && bool(convergence), "Failed writing benchmark data");
    std::ofstream summary(output / "summary.txt");
    summary << std::setprecision(10) << "Passed " << checks << " analytical and diagnostic checks.\n"
            << "Maximum relative error: " << maximumError << "\n"
            << "Convergence error ratio on halving h: " << errorRatio << "\n";
    check(bool(summary), "Cannot write benchmark summary");
    std::cout << "Passed " << checks << " checks; maximum relative error " << maximumError
            << "; convergence ratio " << errorRatio << '\n';
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
