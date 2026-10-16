#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <raytracing/raytracing.hpp>

namespace py = pybind11;
using namespace py::literals;

void bind_common(py::module_ &m);
void bind_point_ray(py::module_ &m);
void bind_beam(py::module_ &m);
