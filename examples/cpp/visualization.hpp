#pragma once
#include <RayTracing/Solver.hpp>
#include <filesystem>

// Application adapters only; neither graphics library enters the solver API.
void plot_impulse(const ImpulseResponse &response, const std::filesystem::path &file);
void view_rays(const rt::SimulationResult &result, bool interactive,
               const std::filesystem::path &screenshot);
