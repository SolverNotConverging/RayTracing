#include <RayTracing/Solver.hpp>
#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

bool hits_bottom(const RefinementResult &path) {
    return std::any_of(path.geometry_.reflections_.begin(), path.geometry_.reflections_.end(),
                      [](const auto &hit) { return hit.surfaceIndex_ == 1; });
}

void check_geometry(const rt::SimulationResult &result) {
    for (const auto &ray : result.rays) {
        const auto &path = ray.refinement;
        check(path.launchDirection_.z() < -0.1, "Tube return did not launch downward");
        check(hits_bottom(path), "Open-rim reflection entered the solved returns");
        check(path.geometry_.receiver_->missDistance_ <= result.settings.refinement.receiverTolerance_,
              "Return missed the point receiver");
        for (const auto &hit : path.geometry_.reflections_) {
            check(hit.position_.z() >= -0.5 - 1e-10 && hit.position_.z() < -1e-5,
                  "Reflection left the physical tube");
            if (hit.surfaceIndex_ == 0)
                check(std::abs(hit.position_.head<2>().norm() - 0.127) < 1e-10,
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
        earlyPass = earlyPass || candidate.reflections_.size() == 1;
        if (candidate.reflections_.size() != 9) return;
        std::vector<std::size_t> sequence;
        for (const auto &hit : candidate.reflections_) sequence.push_back(hit.surfaceIndex_);
        const auto refined = refine_path(origin, launch, origin, surfaces, sequence);
        check(refined.status_ == RefinementStatus::Converged && hits_bottom(refined),
              "Later wall/bottom return was not refined");
        check(std::abs(refined.geometry_.receiver_->pathDistance_ -
                       std::hypot(1.0, 16 * 0.127)) < 1e-10,
              "Unfolded tube return length is incorrect");
        fullReturn = true;
    });
    check(earlyPass && fullReturn, "Reception sphere truncated a later physical return");

    const auto rim = refine_path(origin, Vec3::UnitX(), origin, surfaces, {0});
    check(rim.status_ == RefinementStatus::UnresolvedEdge,
          "Open cylinder rim treated as a smooth specular surface");
    check(deduplicate_paths({rim}).empty(), "Rim reflection entered verified paths");

    const auto patch = rt::load_farfield(TEST_PATTERN);
    rt::SolverConfig config;
    config.rayCount = 1000;
    const auto centered = rt::solve(tube, {origin, patch}, {origin, patch}, config);
    check_geometry(centered);
    check(centered.rays.size() > 1, "Centered wall/bottom return families were lost");
    for (const auto &solved : centered.rays) {
        if (solved.refinement.geometry_.reflections_.size() > 1)
            check(solved.spreading.status_ == SpreadingStatus::Caustic && !solved.field,
                  "On-axis cylindrical focus produced a spurious finite field");
    }
    check(centered.responseRayIndices.size() == 1, "Centered caustic polluted the impulse response");

    const Vec3 shifted(0.005, 0, 0);
    const auto offAxis = rt::solve(tube, {shifted, patch}, {shifted, patch}, config);
    check_geometry(offAxis);
    check(offAxis.rays.size() == 13 && offAxis.responseRayIndices.size() == 13,
          "Shifted tube did not resolve the 13 distinct finite-field returns");
    double peak = 0;
    for (const auto &solved : offAxis.rays) {
        check(solved.field && solved.spreading.status_ == SpreadingStatus::Valid,
              "Shifted return still lies at a receiver caustic");
        peak = std::max(peak, solved.field->receiverField_.norm());
    }
    check(peak > 0 && peak < 300, "Shifted tube has a divergent field");

    // Doubling launches must discover the same discrete roots, not add copies to h(t).
    config.rayCount *= 2;
    const auto denser = rt::solve(tube, {shifted, patch}, {shifted, patch}, config);
    check_geometry(denser);
    check(denser.rays.size() == offAxis.rays.size(), "Refinement left duplicate focused paths");
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
