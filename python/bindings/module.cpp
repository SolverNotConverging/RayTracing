#include "bindings.hpp"

PYBIND11_MODULE(_core, m) {
    m.doc() = "Electromagnetic ray tracing and numerical result persistence.";
    m.attr("__version__") = RAYTRACING_VERSION;
    m.attr("SPEED_OF_LIGHT") = C;
    bind_geometry(m);
    bind_antennas(m);
    bind_results(m);
    bind_solver(m);
}
