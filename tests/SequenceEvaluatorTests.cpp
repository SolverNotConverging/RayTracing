#include "SequenceEvaluator.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double actual, double expected, const char *message) {
    check(std::abs(actual - expected) < 1e-10, message);
}
void near(const Vec3 &actual, const Vec3 &expected, const char *message) {
    check((actual - expected).norm() < 1e-10, message);
}
template<class Action> void rejects(Action action) {
    bool rejected = false;
    try { action(); } catch (const std::invalid_argument &) { rejected = true; }
    check(rejected, "Expected invalid arguments to be rejected");
}
}

int main() try {
    const Vec3 O = Vec3::Zero(), X = Vec3::UnitX(), Y = Vec3::UnitY(), Z = Vec3::UnitZ();
    const std::vector<Surface> wall{Rectangle{Vec3(1, 0, 0), Y, Z, 3, 3}};
    for (const auto &surface : wall) validate_surface(surface);

    const auto direct = evaluate_sequence(O, X, Vec3(3, 0, 0), {}, {});
    check(direct.status_ == SequenceStatus::Valid && direct.receiver_, "Direct ray failed");
    check(direct.reflections_.empty(), "Direct ray acquired reflections");
    near(direct.receiver_->missDistance_, 0, "Direct miss distance incorrect");
    near(direct.receiver_->pathDistance_, 3, "Direct length incorrect");
    const auto miss = evaluate_sequence(O, X, Vec3(3, 4, 0), {}, {});
    check(miss.status_ == SequenceStatus::Valid, "Off-axis valid trial rejected");
    near(miss.receiver_->closestPoint_, Vec3(3, 0, 0), "Closest point incorrect");
    near(miss.receiver_->residual_, Vec3(0, -4, 0), "Residual sign incorrect");
    near(miss.receiver_->missDistance_, 4, "Miss distance incorrect");
    near(miss.receiver_->pathDistance_, 3, "An artificial segment to Rx was added");
    const auto behind = evaluate_sequence(O, X, -X, {}, {});
    check(behind.status_ == SequenceStatus::ReceiverBehind && !behind.receiver_, "Behind receiver accepted");
    const auto coincident = evaluate_sequence(O, X, O, {}, {});
    check(coincident.status_ == SequenceStatus::Valid, "Zero final segment rejected");
    near(coincident.receiver_->pathDistance_, 0, "Zero final segment has length");

    // Tx=(0,0), mirror x=1, Rx=(0,2): exact reflection at (1,1).
    const auto exact = evaluate_sequence(O, Vec3(1, 1, 0).normalized(), Vec3(0, 2, 0), wall, {0});
    check(exact.status_ == SequenceStatus::Valid && exact.reflections_.size() == 1, "Single bounce failed");
    near(exact.reflections_[0].position_, Vec3(1, 1, 0), "Reflection point incorrect");
    near(exact.reflections_[0].normal_, X, "Reflection normal incorrect");
    near(exact.reflections_[0].segmentDistance_, std::sqrt(2.0), "Reflection segment incorrect");
    near(exact.finalDirection_, Vec3(-1, 1, 0).normalized(), "Reflected direction incorrect");
    near(exact.receiver_->missDistance_, 0, "Exact reflection missed Rx");
    near(exact.receiver_->pathDistance_, 2 * std::sqrt(2.0), "Reflected path length incorrect");
    const auto perturbed = evaluate_sequence(O, Vec3(1, 0.9, 0).normalized(), Vec3(0, 2, 0), wall, {0});
    check(perturbed.status_ == SequenceStatus::Valid && perturbed.receiver_->missDistance_ > 0.01,
          "Perturbed direction lost its usable residual");
    near(perturbed.receiver_->residual_.dot(perturbed.finalDirection_), 0, "Residual not perpendicular to final ray");
    const auto missed = evaluate_sequence(O, Y, Vec3(0, 2, 0), wall, {0});
    check(missed.status_ == SequenceStatus::MissedSurface && !missed.receiver_, "Missing reflection not reported");
    const auto pastWall = evaluate_sequence(O, X, Vec3(2, 0, 0), wall, {0});
    check(pastWall.status_ == SequenceStatus::ReceiverBehind, "Receiver behind reflected ray accepted");

    auto obstructed = wall;
    obstructed.push_back(Disk{Vec3(0.5, 0, 0), X, 3});
    for (const auto &surface : obstructed) validate_surface(surface);
    const auto wrong = evaluate_sequence(O, X, -X, obstructed, {0});
    check(wrong.status_ == SequenceStatus::UnexpectedSurface && wrong.blockingSurfaceIndex_ == 1,
          "Unexpected nearer surface ignored");
    check(wrong.reflections_.empty(), "Unexpected surface recorded as requested reflection");
    const auto directBlocked = evaluate_sequence(O, X, Vec3(2, 0, 0), wall, {});
    check(directBlocked.status_ == SequenceStatus::FinalSegmentBlocked && directBlocked.blockingSurfaceIndex_ == 0,
          "Direct obstruction ignored");
    check(directBlocked.receiver_.has_value(), "Blocked approach diagnostics missing");
    auto finalBlockedScene = wall;
    finalBlockedScene.push_back(Sphere{Vec3(0.5, 1.5, 0), 0.1});
    for (const auto &surface : finalBlockedScene) validate_surface(surface);
    const auto finalBlocked = evaluate_sequence(O, Vec3(1, 1, 0).normalized(), Vec3(0, 2, 0), finalBlockedScene, {0});
    check(finalBlocked.status_ == SequenceStatus::FinalSegmentBlocked && finalBlocked.blockingSurfaceIndex_ == 1,
          "Obstruction after reflection ignored");
    check(finalBlocked.reflections_.size() == 1, "Successful prefix lost");
    finalBlockedScene[1] = Sphere{Vec3(-0.5, 2.5, 0), 0.1};
    validate_surface(finalBlockedScene[1]);
    const auto beyond = evaluate_sequence(O, Vec3(1, 1, 0).normalized(), Vec3(0, 2, 0), finalBlockedScene, {0});
    check(beyond.status_ == SequenceStatus::Valid, "Surface beyond receiver incorrectly blocks path");
    const auto endpointBlocked = evaluate_sequence(O, X, X, wall, {});
    check(endpointBlocked.status_ == SequenceStatus::FinalSegmentBlocked, "Surface at endpoint should count as blocking");

    const std::vector<Surface> corridor{wall[0], Rectangle{Vec3(-1, 0, 0), Y, Z, 3, 3}};
    for (const auto &surface : corridor) validate_surface(surface);
    const auto repeated = evaluate_sequence(O, X, O, corridor, {0, 1, 0});
    check(repeated.status_ == SequenceStatus::Valid && repeated.reflections_.size() == 3,
          "Repeated surface sequence failed");
    near(repeated.receiver_->pathDistance_, 6, "Multiple reflection distance incorrect");
    const auto prefixFailure = evaluate_sequence(O, X, O, corridor, {0, 0});
    check(prefixFailure.status_ == SequenceStatus::UnexpectedSurface && prefixFailure.reflections_.size() == 1,
          "Failure did not retain replayed prefix");

    const std::vector<Surface> curved{Sphere{Vec3(2, 0, 0), 1}};
    validate_surface(curved[0]);
    const auto curvedHit = evaluate_sequence(O, X, -X, curved, {0});
    check(curvedHit.status_ == SequenceStatus::Valid, "Curved sequence failed");
    near(curvedHit.reflections_[0].normal_, -X, "Curved normal incorrect");
    near(curvedHit.receiver_->pathDistance_, 3, "Curved path length incorrect");

    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(0.4, Vec3(1, 2, 3).normalized()).toRotationMatrix();
    const Vec3 shift(3, -2, 1);
    const std::vector<Surface> rotated{Rectangle{rotation * X + shift, rotation * Y, rotation * Z, 3, 3}};
    validate_surface(rotated[0]);
    const auto transformed = evaluate_sequence(shift, rotation * Vec3(1, 1, 0).normalized(),
        rotation * Vec3(0, 2, 0) + shift, rotated, {0});
    check(transformed.status_ == SequenceStatus::Valid, "Rotated sequence failed");
    near(transformed.receiver_->pathDistance_, exact.receiver_->pathDistance_, "Rigid transform changed length");
    near(transformed.receiver_->missDistance_, 0, "Rigid transform introduced miss");

    const Surface cornerX = Rectangle{X, Y, Z, 1, 1};
    const Surface cornerY = Rectangle{Y, X, Z, 1, 1};
    validate_surface(cornerX);
    validate_surface(cornerY);
    for (const auto &corner : {std::vector<Surface>{cornerX, cornerY}, std::vector<Surface>{cornerY, cornerX}}) {
        const auto ambiguous = evaluate_sequence(O, Vec3(1, 1, 0).normalized(), O, corner, {0});
        check(ambiguous.status_ == SequenceStatus::AmbiguousHit && ambiguous.reflections_.empty(),
              "Sequence accepted arbitrary corner normal");
        near(ambiguous.finalOrigin_, Vec3(1, 1, 0), "Corner position not saved");
        const auto blockedCorner = evaluate_sequence(O, Vec3(1, 1, 0).normalized(), Vec3(2, 2, 0), corner, {});
        check(blockedCorner.status_ == SequenceStatus::AmbiguousHit, "Final-leg ambiguous blocker not reported");
    }
    const std::vector<Surface> seam{
        Triangle{Vec3(1, -1, -1), Vec3(1, 1, -1), Vec3(1, 1, 1)},
        Triangle{Vec3(1, -1, -1), Vec3(1, -1, 1), Vec3(1, 1, 1)}};
    for (const auto &surface : seam) validate_surface(surface);
    for (std::size_t index : {0u, 1u}) {
        const auto seamResult = evaluate_sequence(O, X, -X, seam, {index});
        check(seamResult.status_ == SequenceStatus::Valid && seamResult.reflections_[0].surfaceIndex_ == index,
              "Equivalent seam face rejected by sequence evaluator");
    }

    rejects([&] { evaluate_sequence(O, X, X, wall, {1}); });
    rejects([&] { evaluate_sequence(O, Vec3::Zero(), X, wall, {}); });
    rejects([&] { evaluate_sequence(O, 2 * X, X, wall, {}); });
    rejects([&] { evaluate_sequence(O, X, X, wall, {}, -1); });
    rejects([&] { evaluate_sequence(O, X, Vec3(std::numeric_limits<double>::quiet_NaN(), 0, 0), wall, {}); });
    std::cout << "All sequence evaluator checks passed\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
