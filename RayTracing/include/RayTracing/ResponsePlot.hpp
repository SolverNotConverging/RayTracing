#pragma once
#include "ImpulseResponse.hpp"
#include <filesystem>

// Legacy wrapper exporting a magnitude/phase PNG. Use rt::save_impulse_csv for data.
// Optional interactive display occurs immediately after export.
void plot_impulse_response(const ImpulseResponse &response,
                           const std::vector<ReconstructedField> &fields, const std::filesystem::path &outputDirectory,
                           bool showWindow = true);
