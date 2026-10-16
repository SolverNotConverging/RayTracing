#include "bindings.hpp"

PYBIND11_MODULE(_core, m) {
    m.doc() = "Native core of the raytracing package";
    m.attr("__version__") = RAYTRACING_VERSION;
    bind_common(m);
    bind_point_ray(m);
    bind_beam(m);
}
