#include "TraceConfig.hpp"

#include <fstream>
#include <limits>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

TraceOptions load_trace_options(const std::filesystem::path &path) {
    try {
        std::ifstream input(path);
        if (!input) {
            throw std::runtime_error("Cannot open file");
        }

        const json config = json::parse(input);
        if (!config.is_object()) {
            throw std::runtime_error("Expected a JSON object");
        }

        // at() reports a missing key; check types before converting values.
        const auto &reflections = config.at("maxReflections");
        const auto &distance = config.at("maxDistance");
        const auto &index = config.at("refractiveIndex");
        if (!reflections.is_number_integer() ||
            reflections.get<double>() < 0.0 ||
            reflections.get<double>() > std::numeric_limits<int>::max()) {
            throw std::runtime_error("maxReflections must be a nonnegative integer within the int range");
        }
        if (!distance.is_number() || !index.is_number()) {
            throw std::runtime_error("maxDistance and refractiveIndex must be numbers");
        }

        TraceOptions options;
        options.maxReflections = reflections.get<int>();
        options.maxDistance = distance.get<double>();
        options.refractiveIndex = index.get<double>();
        if (!std::isfinite(options.maxDistance) || options.maxDistance <= 0.0) {
            throw std::runtime_error("maxDistance must be finite and positive (metres)");
        }
        if (!std::isfinite(options.refractiveIndex) || options.refractiveIndex <= 0.0) {
            throw std::runtime_error("refractiveIndex must be finite and positive");
        }
        return options;
    } catch (const std::exception &error) {
        throw std::runtime_error("Cannot load trace options from '" + path.string() + "': " + error.what());
    }
}
