#include <raytracing/refinement.hpp>
#include <raytracing/tracing.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace rt;

namespace {
void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
void converged(const RefinementResult &result) {
    if (result.status != RefinementStatus::Converged)
        throw std::runtime_error("Refinement status " + std::to_string(static_cast<int>(result.status)));
    check(result.geometry.status == SequenceStatus::Valid && result.geometry.receiver &&
          result.geometry.receiver->missDistance <= 1e-6, "Converged path missed Rx");
    check(std::abs(result.launchDirection.norm() - 1) < 1e-12, "Direction lost normalization");
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
    check((direct.launchDirection - X).norm() < 1e-6, "Direct solution incorrect");

    const auto reflected = refine_path(O, seed, rx, wall, {0});
    converged(reflected);
    check((reflected.geometry.reflections[0].position - Vec3(1, 1, 0)).norm() < 1e-6,
          "Single reflection did not reach analytical solution");
    check(std::abs(reflected.geometry.receiver->pathDistance - 2 * std::sqrt(2.0)) < 1e-6,
          "Reflected path length incorrect");
    check(reflected.iterations > 0, "Off-axis seed was not refined");

    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(0.7, Vec3(1, 2, 3).normalized()).toRotationMatrix();
    const Vec3 shift(3, -2, 1);
    const std::vector<Surface> rotated{Rectangle{shift + rotation * X, rotation * Y, rotation * Z, 3, 3}};
    validate_surface(rotated[0]);
    const auto transformed = refine_path(shift, rotation * seed, shift + rotation * rx, rotated, {0});
    converged(transformed);
    check((transformed.geometry.reflections[0].position - (shift + rotation * Vec3(1, 1, 0))).norm() < 1e-6,
          "Refinement depends on coordinate orientation");

    const std::vector<Surface> corridor{wall[0], Rectangle{-X, Y, Z, 3, 3}};
    for (const auto &surface : corridor) validate_surface(surface);
    const auto multiple = refine_path(O, Vec3(1, 0.02, -0.01).normalized(), O, corridor, {0, 1, 0});
    converged(multiple);
    check(std::abs(multiple.geometry.receiver->pathDistance - 6) < 1e-6,
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
        check((curved.geometry.reflections[0].position - hit).norm() < 1e-6,
              "Curved solution incorrect");
        check(std::abs(curved.geometry.receiver->pathDistance - 5) < 1e-6,
              "Curved length incorrect");
    }

    RefinementOptions noIterations;
    noIterations.maxIterations = 0;
    const auto unfinished = refine_path(O, seed, rx, wall, {0}, noIterations);
    check(unfinished.status == RefinementStatus::IterationLimit, "Iteration limit ignored");
    converged(refine_path(O, X, 3 * X, {}, {}, noIterations));
    RefinementOptions shortPath;
    shortPath.maxPathDistance = 2;
    const auto tooLong = refine_path(O, seed, rx, wall, {0}, shortPath);
    check(tooLong.status == RefinementStatus::PathTooLong, "Distance limit ignored");
    check(refine_path(O, Y, rx, wall, {0}).status == RefinementStatus::InvalidGeometry,
          "Missing required surface accepted");
    check(refine_path(O, X, -X, {}, {}).status == RefinementStatus::InvalidGeometry,
          "Backward receiver accepted");
    check(refine_path(O, X, 2 * X, wall, {}).status == RefinementStatus::InvalidGeometry,
          "Direct obstruction ignored");
    auto blocked = wall;
    blocked.push_back(Sphere{Vec3(0.5, 1.5, 0), 0.1});
    validate_surface(blocked.back());
    check(refine_path(O, Vec3(1, 1, 0).normalized(), rx, blocked, {0}).status == RefinementStatus::InvalidGeometry,
          "Final reflected segment obstruction ignored");
    // The unconstrained reflection point is outside this finite rectangle.
    const std::vector<Surface> smallWall{Rectangle{X, Y, Z, 0.5, 3}};
    validate_surface(smallWall[0]);
    check(refine_path(O, Vec3(1, 0.4, 0).normalized(), rx, smallWall, {0}).status != RefinementStatus::Converged,
          "Optimizer bypassed finite surface boundary");

    // Same physical path reached from separate seeds; best residual is retained.
    const auto second = refine_path(O, Vec3(1, 1.05, -0.03).normalized(), rx, wall, {0});
    converged(second);
    const std::vector<RefinementResult> duplicates{reflected, second, unfinished, tooLong};
    const auto unique = deduplicate_paths(duplicates);
    check(unique.size() == 1, "Duplicates or failed paths retained");
    check(duplicates[unique[0]].geometry.receiver->missDistance <=
          duplicates[1 - unique[0]].geometry.receiver->missDistance, "Best representative not retained");
    const auto opposite = refine_path(O, Vec3(-1, 0.9, 0.02).normalized(), rx, corridor, {1});
    converged(opposite);
    check(deduplicate_paths({reflected, opposite}).size() == 2, "Distinct reflection sequences merged");

    // Comparator fixtures: sharing a surface sequence alone does not imply duplication.
    auto differentPoint = reflected;
    differentPoint.geometry.reflections[0].position.z() += 1e-3;
    auto differentArrival = reflected;
    differentArrival.geometry.finalDirection = (reflected.geometry.finalDirection + 1e-3 * Z).normalized();
    auto differentLength = reflected;
    differentLength.geometry.receiver->pathDistance += 1e-3;
    check(deduplicate_paths({reflected, differentPoint, differentArrival, differentLength}).size() == 4,
          "Different solutions with a common sequence merged");
    auto a = reflected, b = reflected, c = reflected;
    a.geometry.receiver->missDistance = 0;
    b.geometry.receiver->missDistance = 1e-9;
    c.geometry.receiver->missDistance = 2e-9;
    b.geometry.reflections[0].position.z() += 0.75e-5;
    c.geometry.reflections[0].position.z() += 1.5e-5;
    check(deduplicate_paths({a, b, c}).size() == 2, "Transitive duplicate merging swallowed a distinct path");

    // Recorded events refer to actual completed reflections and incident samples.
    TracedRay ray(O, Vec3(1, 1, 0).normalized());
    TraceOptions traceOptions;
    const auto traced = trace_ray(ray, wall, rx, 0.1, traceOptions);
    check(traced.status == TraceStatus::Received && ray.reflections.size() == 1, "Reflection event missing");
    const auto &event = ray.reflections[0];
    check(event.surfaceIndex == 0 && (event.normal - X).norm() < 1e-12, "Event surface or normal incorrect");
    check(event.vertexIndex < ray.vertices.size() &&
          (ray.vertices[event.vertexIndex] - Vec3(1, 1, 0)).norm() < 1e-12,
          "Event hit vertex incorrect");
    TracedRay limited(O, X);
    traceOptions.maxReflections = 0;
    trace_ray(limited, wall, -X, 0.1, traceOptions);
    check(limited.reflections.empty(), "Unperformed reflection recorded");

    const std::vector<Surface> corner{Rectangle{X, Y, Z, 1, 1}, Rectangle{Y, X, Z, 1, 1}};
    for (const auto &surface : corner) validate_surface(surface);
    const auto ambiguous = refine_path(O, Vec3(1, 1, 0).normalized(), O, corner, {0});
    check(ambiguous.status == RefinementStatus::UnresolvedCorner &&
          ambiguous.geometry.status == SequenceStatus::AmbiguousHit, "Refinement accepted corner");
    check(deduplicate_paths({ambiguous}).empty(), "Ambiguous candidate retained");

    // Monostatic reception: near-corner rays leave the Rx sphere, reflect twice,
    // and return into it. Refinement collapses the two hits toward the corner.
    for (double offset : {-0.01, 0.01}) {
        const Vec3 launch = Vec3(1, 1 + offset, 0).normalized();
        TracedRay candidate(O, launch);
        const auto received = trace_ray(candidate, corner, O, 0.1);
        check(received.status == TraceStatus::Received && received.reflections == 2 && received.distance > 2,
              "Monostatic corner return terminated at launch or failed reception");
        check(std::abs(candidate.position().norm() - 0.1) < 1e-10,
              "Return did not stop at sphere entry");
        std::vector<std::size_t> sequence;
        for (const auto &event : candidate.reflections) sequence.push_back(event.surfaceIndex);
        const auto refined = refine_path(O, launch, O, corner, sequence);
        check(refined.status == RefinementStatus::UnresolvedCorner,
              "Collapsing corner return accepted as an ordinary point-Rx solution");
        check(deduplicate_paths({refined}).empty(), "Unresolved corner entered verified solutions");
        check(candidate.reflections.size() == 2, "Refinement destroyed coarse return");
    }
    // A nearby bistatic receiver has a genuine pair of separated reflection points.
    const auto separated = refine_path(O, Vec3(1, 0.89, 0).normalized(), Vec3(0, 0.2, 0), corner, {0, 1});
    converged(separated);
    check(separated.geometry.reflections[1].segmentDistance > 0.1,
          "Separated corner solution unexpectedly collapsed");
    TracedRay freeSpace(O, X);
    check(trace_ray(freeSpace, {}, O, 0.1).status == TraceStatus::Escaped,
          "Launch inside Rx without a return was received");
    TracedRay normalReturn(O, X);
    const auto singleReturn = trace_ray(normalReturn, wall, O, 0.1);
    check(singleReturn.status == TraceStatus::Received && singleReturn.reflections == 1,
          "Single-face monostatic return was lost");

    std::cout << "All path refinement checks passed\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
