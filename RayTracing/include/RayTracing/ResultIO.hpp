#pragma once
#include "Solver.hpp"

namespace rt {
    void save_h5(const SimulationResult &result, const std::filesystem::path &path);

    SimulationResult load_h5(const std::filesystem::path &path);

    void save_impulse_csv(const ImpulseResponse &response,
                          const std::filesystem::path &path);

    ImpulseResponse load_impulse_csv(const std::filesystem::path &path);
} // namespace rt
