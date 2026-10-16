#include <raytracing/field.hpp>
#include <raytracing/tracing.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rt {
    PathField reconstruct_field(const RefinementResult &path, const SpreadingResult &spreading,
                                const Vec3C &sourceField, const std::vector<Material> &surfaceMaterials,
                                const FieldOptions &options) {
        if (path.status != RefinementStatus::Converged || spreading.status != SpreadingStatus::Valid)
            throw std::invalid_argument("Field reconstruction needs a converged path with valid spreading");
        const auto &geometry = path.geometry;
        const double f = options.frequencyHz;
        PathField result;
        Vec3C field = sourceField;
        Vec3 direction = path.launchDirection;
        double length = geometry.receiver->finalSegmentDistance;
        for (const auto &hit: geometry.reflections) {
            const double c = std::min(1.0, std::abs(direction.dot(hit.normal)));
            const auto r = fresnel(options.background, surfaceMaterials.at(hit.surfaceIndex), c, f);
            field = reflection_matrix(direction, hit.normal, r.rs, r.rp) * field;
            result.reflections.push_back({std::acos(c), r.rs, r.rp});
            direction = reflected_direction(direction, hit.normal).normalized();
            length += hit.segmentDistance;
        }

        result.pathLength = length;
        result.opticalPath = options.background.index(f).real() * length;
        result.delay = result.opticalPath / speedOfLight;
        result.frequencyHz = f;
        result.fieldFactor = *spreading.fieldFactor;
        result.causticCount = spreading.causticCount;
        result.arrivalDirection = direction;
        result.transportedField = field;
        // Reduce the cycle count before forming the phase to keep it accurate for long paths.
        const Complex phase = std::polar(1.0, 2.0 * pi * std::remainder(f * result.delay, 1.0));
        const Complex maslov[] = {{1, 0}, {0, -1}, {-1, 0}, {0, 1}};
        result.field = result.fieldFactor * phase * maslov[result.causticCount % 4] * field;
        result.normalizedField = result.field / options.sourceReference;
        return result;
    }

    ImpulseResponse group_taps(std::vector<ImpulseTap> paths, double frequencyHz, double delayTolerance) {
        std::stable_sort(paths.begin(), paths.end(), [](const auto &a, const auto &b) { return a.delay < b.delay; });
        ImpulseResponse response{frequencyHz, {}};
        for (auto &path: paths) {
            if (response.taps.empty() || path.delay - response.taps.back().delay > delayTolerance) {
                response.taps.push_back(std::move(path));
                continue;
            }
            auto &tap = response.taps.back();
            tap.coefficient += path.coefficient;
            tap.field += path.field;
            tap.paths.insert(tap.paths.end(), path.paths.begin(), path.paths.end());
        }
        return response;
    }

    ImpulseResponse calculate_impulse_response(const std::vector<PathField> &fields, const Vec3C &analyzer,
                                               double delayTolerance) {
        std::vector<ImpulseTap> paths;
        for (std::size_t i = 0; i < fields.size(); ++i)
            paths.push_back({fields[i].delay, analyzer.dot(fields[i].normalizedField), fields[i].normalizedField, {i}});
        return group_taps(std::move(paths), fields.empty() ? 0.0 : fields.front().frequencyHz, delayTolerance);
    }

    Complex frequency_response(const ImpulseResponse &response, double offsetHz) {
        Complex sum = 0.0;
        for (const auto &tap: response.taps)
            sum += tap.coefficient * std::polar(1.0, 2.0 * pi * std::remainder(offsetHz * tap.delay, 1.0));
        return sum;
    }
} // namespace rt
