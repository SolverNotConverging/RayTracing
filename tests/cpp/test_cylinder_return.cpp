#include <raytracing/solver.hpp>
#include <algorithm>
#include <iostream>
#include <stdexcept>

using namespace rt;

namespace {
void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

bool hits_bottom(const RefinementResult &path) {
    return std::any_of(path.geometry.reflections.begin(), path.geometry.reflections.end(),
                      [](const auto &hit) { return hit.surfaceIndex == 1; });
}

void check_geometry(const rt::SimulationResult &result) {
    for (const auto &ray : result.paths) {
        const auto &path = ray.refinement;
        check(path.launchDirection.z() < -0.1, "Tube return did not launch downward");
        check(hits_bottom(path), "Open-rim reflection entered the solved returns");
        check(path.geometry.receiver->missDistance <= result.config.refinement.receiverTolerance,
              "Return missed the point receiver");
        for (const auto &hit : path.geometry.reflections) {
            check(hit.position.z() >= -0.5 - 1e-10 && hit.position.z() < -1e-5,
                  "Reflection left the physical tube");
            if (hit.surfaceIndex == 0)
                check(std::abs(hit.position.head<2>().norm() - 0.127) < 1e-10,
                      "Wall reflection left the cylinder");
        }
    }
}
}

int main() try {
    const Vec3 origin = Vec3::Zero();
    const auto pec = rt::Material::pec();
    rt::Scene tube;
    tube.add(Cylinder{{0, 0, -0.25}, Vec3::UnitZ(), 0.127, 0.25, false}, pec);
    tube.add(Disk{{0, 0, -0.5}, -Vec3::UnitZ(), 0.127}, pec);
    const auto &surfaces = tube.surfaces();

    // An exact eight-wall + bottom return passes near Rx before reaching bottom.
    // The capture sphere must not absorb it at that first, unrefinable near-pass.
    const Vec3 launch = Vec3(8 * 0.127, 0, -0.5).normalized();
    TracedRay ray(origin, launch);
    TraceOptions tracing{9, 20};
    bool earlyPass = false, fullReturn = false;
    trace_ray(ray, surfaces, origin, 0.12, tracing, [&](const TracedRay &candidate) {
        earlyPass = earlyPass || candidate.reflections.size() == 1;
        if (candidate.reflections.size() != 9) return;
        std::vector<std::size_t> sequence;
        for (const auto &hit : candidate.reflections) sequence.push_back(hit.surfaceIndex);
        const auto refined = refine_path(origin, launch, origin, surfaces, sequence);
        check(refined.status == RefinementStatus::Converged && hits_bottom(refined),
              "Later wall/bottom return was not refined");
        check(std::abs(refined.geometry.receiver->pathDistance -
                       std::hypot(1.0, 16 * 0.127)) < 1e-10,
              "Unfolded tube return length is incorrect");
        fullReturn = true;
    });
    check(earlyPass && fullReturn, "Reception sphere truncated a later physical return");

    const auto rim = refine_path(origin, Vec3::UnitX(), origin, surfaces, {0});
    check(rim.status == RefinementStatus::UnresolvedEdge,
          "Open cylinder rim treated as a smooth specular surface");
    check(deduplicate_paths({rim}).empty(), "Rim reflection entered verified paths");

    const auto patch = rt::load_farfield(TEST_DATA "/simple_patch.ffs");
    rt::SolverConfig config;
    config.rayCount = 1000;
    const auto centered = rt::solve(tube, {origin, patch}, {origin, patch}, config);
    check_geometry(centered);
    check(centered.paths.size() > 1, "Centered wall/bottom return families were lost");
    for (const auto &solved : centered.paths) {
        if (solved.refinement.geometry.reflections.size() > 1)
            check(solved.spreading.status == SpreadingStatus::Caustic && !solved.field,
                  "On-axis cylindrical focus produced a spurious finite field");
    }
    check(centered.impulseResponse.taps.size() == 1, "Centered caustic polluted the impulse response");

    const Vec3 shifted(0.005, 0, 0);
    const auto offAxis = rt::solve(tube, {shifted, patch}, {shifted, patch}, config);
    check_geometry(offAxis);
    check(offAxis.paths.size() == 13,
          "Shifted tube did not resolve the 13 distinct finite-field returns");
    double peak = 0;
    for (const auto &solved : offAxis.paths) {
        check(solved.field && solved.spreading.status == SpreadingStatus::Valid,
              "Shifted return still lies at a receiver caustic");
        peak = std::max(peak, solved.field->field.norm());
    }
    check(peak > 0 && peak < 300, "Shifted tube has a divergent field");
    bool internalFocus = false;
    for (const auto &solved : offAxis.paths)
        internalFocus = internalFocus || solved.field->causticCount > 0;
    check(internalFocus, "Shifted tube did not exercise an internal focus");

    // Doubling launches must discover the same discrete roots, not add copies to h(t).
    config.rayCount *= 2;
    const auto denser = rt::solve(tube, {shifted, patch}, {shifted, patch}, config);
    check_geometry(denser);
    check(denser.paths.size() == offAxis.paths.size(), "Refinement left duplicate focused paths");
    const Complex response = frequency_response(offAxis.impulseResponse);
    check(std::abs(frequency_response(denser.impulseResponse) - response) <
              1e-5 * std::max(1.0, std::abs(response)),
          "Shifted response depends on the number of launch seeds");
    std::cout << "Cylinder returns, caustic rejection, shifted fields and sampling convergence passed\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
