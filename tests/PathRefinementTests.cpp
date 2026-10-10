#include "PathRefinement.hpp"
#include "Tracing.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
void converged(const RefinementResult &result) {
    if (result.status_ != RefinementStatus::Converged)
        throw std::runtime_error(refinement_status_name(result.status_));
    check(result.geometry_.status_ == SequenceStatus::Valid && result.geometry_.receiver_ &&
          result.geometry_.receiver_->missDistance_ <= 1e-6, "Converged path missed Rx");
    check(std::abs(result.launchDirection_.norm() - 1) < 1e-12, "Direction lost normalization");
}
template<class Action> void rejects(Action action) {
    bool rejected = false;
    try { action(); } catch (const std::invalid_argument &) { rejected = true; }
    check(rejected, "Invalid arguments accepted");
}
}

int main() try {
    const Vec3 O = Vec3::Zero(), X = Vec3::UnitX(), Y = Vec3::UnitY(), Z = Vec3::UnitZ();
    const Vec3 rx(0, 2, 0), seed = Vec3(1, 0.9, 0.02).normalized();
    const std::vector<Surface> wall{Rectangle{X, Y, Z, 3, 3}};
    for (const auto &surface : wall) validate_surface(surface);
    const auto direct = refine_path(O, Vec3(1, 0.1, -0.1).normalized(), 3 * X, {}, {});
    converged(direct);
    check((direct.launchDirection_ - X).norm() < 1e-6, "Direct solution incorrect");

    const auto reflected = refine_path(O, seed, rx, wall, {0});
    converged(reflected);
    check((reflected.geometry_.reflections_[0].position_ - Vec3(1, 1, 0)).norm() < 1e-6,
          "Single reflection did not reach analytical solution");
    check(std::abs(reflected.geometry_.receiver_->pathDistance_ - 2 * std::sqrt(2.0)) < 1e-6,
          "Reflected path length incorrect");
    check(reflected.iterations_ > 0, "Off-axis seed was not refined");

    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(0.7, Vec3(1, 2, 3).normalized()).toRotationMatrix();
    const Vec3 shift(3, -2, 1);
    const std::vector<Surface> rotated{Rectangle{shift + rotation * X, rotation * Y, rotation * Z, 3, 3}};
    validate_surface(rotated[0]);
    const auto transformed = refine_path(shift, rotation * seed, shift + rotation * rx, rotated, {0});
    converged(transformed);
    check((transformed.geometry_.reflections_[0].position_ - (shift + rotation * Vec3(1, 1, 0))).norm() < 1e-6,
          "Refinement depends on coordinate orientation");

    const std::vector<Surface> corridor{wall[0], Rectangle{-X, Y, Z, 3, 3}};
    for (const auto &surface : corridor) validate_surface(surface);
    const auto multiple = refine_path(O, Vec3(1, 0.02, -0.01).normalized(), O, corridor, {0, 1, 0});
    converged(multiple);
    check(std::abs(multiple.geometry_.receiver_->pathDistance_ - 6) < 1e-6,
          "Repeated reflection solution incorrect");

    // Construct exact curved-surface paths from their known normal at (-1,0,0).
    const Vec3 hit(-1, 0, 0), incident = Vec3(1, 0.3, 0.1).normalized();
    const Vec3 outgoing(-incident.x(), incident.y(), incident.z());
    const Vec3 txCurved = hit - 2 * incident, rxCurved = hit + 3 * outgoing;
    for (const Surface &surface : std::vector<Surface>{Sphere{O, 1}, Cylinder{O, Z, 1, 4, true}}) {
        validate_surface(surface);
        const auto curved = refine_path(txCurved, (incident + Vec3(0, 0.003, -0.002)).normalized(),
                                        rxCurved, {surface}, {0});
        converged(curved);
        check((curved.geometry_.reflections_[0].position_ - hit).norm() < 1e-6,
              "Curved solution incorrect");
        check(std::abs(curved.geometry_.receiver_->pathDistance_ - 5) < 1e-6,
              "Curved length incorrect");
    }

    RefinementOptions noIterations;
    noIterations.maxIterations_ = 0;
    const auto unfinished = refine_path(O, seed, rx, wall, {0}, noIterations);
    check(unfinished.status_ == RefinementStatus::IterationLimit, "Iteration limit ignored");
    converged(refine_path(O, X, 3 * X, {}, {}, noIterations));
    RefinementOptions shortPath;
    shortPath.maxPathDistance_ = 2;
    const auto tooLong = refine_path(O, seed, rx, wall, {0}, shortPath);
    check(tooLong.status_ == RefinementStatus::PathTooLong, "Distance limit ignored");
    check(refine_path(O, Y, rx, wall, {0}).status_ == RefinementStatus::InvalidGeometry,
          "Missing required surface accepted");
    check(refine_path(O, X, -X, {}, {}).status_ == RefinementStatus::InvalidGeometry,
          "Backward receiver accepted");
    check(refine_path(O, X, 2 * X, wall, {}).status_ == RefinementStatus::InvalidGeometry,
          "Direct obstruction ignored");
    auto blocked = wall;
    blocked.push_back(Sphere{Vec3(0.5, 1.5, 0), 0.1});
    validate_surface(blocked.back());
    check(refine_path(O, Vec3(1, 1, 0).normalized(), rx, blocked, {0}).status_ == RefinementStatus::InvalidGeometry,
          "Final reflected segment obstruction ignored");
    // The unconstrained reflection point is outside this finite rectangle.
    const std::vector<Surface> smallWall{Rectangle{X, Y, Z, 0.5, 3}};
    validate_surface(smallWall[0]);
    check(refine_path(O, Vec3(1, 0.4, 0).normalized(), rx, smallWall, {0}).status_ != RefinementStatus::Converged,
          "Optimizer bypassed finite surface boundary");

    // Same physical path reached from separate seeds; best residual is retained.
    const auto second = refine_path(O, Vec3(1, 1.05, -0.03).normalized(), rx, wall, {0});
    converged(second);
    const std::vector<RefinementResult> duplicates{reflected, second, unfinished, tooLong};
    const auto unique = deduplicate_paths(duplicates);
    check(unique.size() == 1, "Duplicates or failed paths retained");
    check(duplicates[unique[0]].geometry_.receiver_->missDistance_ <=
          duplicates[1 - unique[0]].geometry_.receiver_->missDistance_, "Best representative not retained");
    const auto opposite = refine_path(O, Vec3(-1, 0.9, 0.02).normalized(), rx, corridor, {1});
    converged(opposite);
    check(deduplicate_paths({reflected, opposite}).size() == 2, "Distinct reflection sequences merged");

    // Comparator fixtures: sharing a surface sequence alone does not imply duplication.
    auto differentPoint = reflected;
    differentPoint.geometry_.reflections_[0].position_.z() += 1e-3;
    auto differentArrival = reflected;
    differentArrival.geometry_.finalDirection_ = (reflected.geometry_.finalDirection_ + 1e-3 * Z).normalized();
    auto differentLength = reflected;
    differentLength.geometry_.receiver_->pathDistance_ += 1e-3;
    check(deduplicate_paths({reflected, differentPoint, differentArrival, differentLength}).size() == 4,
          "Different solutions with a common sequence merged");
    auto a = reflected, b = reflected, c = reflected;
    a.geometry_.receiver_->missDistance_ = 0;
    b.geometry_.receiver_->missDistance_ = 1e-9;
    c.geometry_.receiver_->missDistance_ = 2e-9;
    b.geometry_.reflections_[0].position_.z() += 0.75e-5;
    c.geometry_.reflections_[0].position_.z() += 1.5e-5;
    check(deduplicate_paths({a, b, c}).size() == 2, "Transitive duplicate merging swallowed a distinct path");

    // Recorded events refer to actual completed reflections and incident samples.
    TracedRay ray(O, Vec3(1, 1, 0).normalized());
    TraceOptions traceOptions;
    const auto traced = trace_ray(ray, wall, rx, 0.1, traceOptions);
    check(traced.status_ == TraceStatus::Received && ray.reflections_.size() == 1, "Reflection event missing");
    const auto &event = ray.reflections_[0];
    check(event.surfaceIndex_ == 0 && (event.normal_ - X).norm() < 1e-12, "Event surface or normal incorrect");
    check(event.vertexIndex_ < ray.vertices_.size() &&
          (ray.vertices_[event.vertexIndex_] - Vec3(1, 1, 0)).norm() < 1e-12,
          "Event hit vertex incorrect");
    TracedRay limited(O, X);
    traceOptions.maxReflections_ = 0;
    trace_ray(limited, wall, -X, 0.1, traceOptions);
    check(limited.reflections_.empty(), "Unperformed reflection recorded");

    const std::vector<Surface> corner{Rectangle{X, Y, Z, 1, 1}, Rectangle{Y, X, Z, 1, 1}};
    for (const auto &surface : corner) validate_surface(surface);
    const auto ambiguous = refine_path(O, Vec3(1, 1, 0).normalized(), O, corner, {0});
    check(ambiguous.status_ == RefinementStatus::UnresolvedCorner &&
          ambiguous.geometry_.status_ == SequenceStatus::AmbiguousHit, "Refinement accepted corner");
    check(deduplicate_paths({ambiguous}).empty(), "Ambiguous candidate retained");

    // Monostatic reception: near-corner rays leave the Rx sphere, reflect twice,
    // and return into it. Refinement collapses the two hits toward the corner.
    for (double offset : {-0.01, 0.01}) {
        const Vec3 launch = Vec3(1, 1 + offset, 0).normalized();
        TracedRay candidate(O, launch);
        const auto received = trace_ray(candidate, corner, O, 0.1);
        check(received.status_ == TraceStatus::Received && received.reflections_ == 2 && received.distance_ > 2,
              "Monostatic corner return terminated at launch or failed reception");
        check(std::abs(candidate.position().norm() - 0.1) < 1e-10,
              "Return did not stop at sphere entry");
        std::vector<std::size_t> sequence;
        for (const auto &event : candidate.reflections_) sequence.push_back(event.surfaceIndex_);
        const auto refined = refine_path(O, launch, O, corner, sequence);
        check(refined.status_ == RefinementStatus::UnresolvedCorner,
              "Collapsing corner return accepted as an ordinary point-Rx solution");
        check(deduplicate_paths({refined}).empty(), "Unresolved corner entered verified solutions");
        check(candidate.reflections_.size() == 2, "Refinement destroyed coarse return");
    }
    // A nearby bistatic receiver has a genuine pair of separated reflection points.
    const auto separated = refine_path(O, Vec3(1, 0.89, 0).normalized(), Vec3(0, 0.2, 0), corner, {0, 1});
    converged(separated);
    check(separated.geometry_.reflections_[1].segmentDistance_ > 0.1,
          "Separated corner solution unexpectedly collapsed");
    TracedRay freeSpace(O, X);
    check(trace_ray(freeSpace, {}, O, 0.1).status_ == TraceStatus::Escaped,
          "Launch inside Rx without a return was received");
    TracedRay normalReturn(O, X);
    const auto singleReturn = trace_ray(normalReturn, wall, O, 0.1);
    check(singleReturn.status_ == TraceStatus::Received && singleReturn.reflections_ == 1,
          "Single-face monostatic return was lost");

    rejects([&] { refine_path(O, 2 * X, rx, wall, {0}); });
    rejects([&] { refine_path(O, X, rx, wall, {2}); });
    RefinementOptions bad;
    bad.receiverTolerance_ = 0;
    rejects([&] { refine_path(O, X, rx, wall, {0}, bad); });
    DeduplicationOptions badDedup;
    badDedup.angleTolerance_ = -1;
    rejects([&] { deduplicate_paths({}, badDedup); });
    std::cout << "All path refinement checks passed\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
