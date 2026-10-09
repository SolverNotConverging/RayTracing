#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <RayTracing/Solver.hpp>

namespace py = pybind11;

void bind_geometry(py::module_ &m);

void bind_antennas(py::module_ &m);

void bind_solver(py::module_ &m);

void bind_results(py::module_ &m);
