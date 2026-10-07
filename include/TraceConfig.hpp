#pragma once

#include "EMRay.hpp"
#include <filesystem>

// All three TraceOptions fields are required in the JSON file.
TraceOptions load_trace_options(const std::filesystem::path &path);
