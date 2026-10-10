#include "bindings.hpp"

PYBIND11_MODULE(_core, m) {
    m.doc() = "Electromagnetic ray tracing and numerical result persistence.";
    m.attr("__version__") = RAYTRACING_VERSION;
    m.attr("SPEED_OF_LIGHT") = rt::constants::speedOfLight;
    m.attr("VACUUM_PERMEABILITY") = rt::constants::vacuumPermeability;
    m.attr("VACUUM_PERMITTIVITY") = rt::constants::vacuumPermittivity;
    m.attr("FREE_SPACE_IMPEDANCE") = rt::constants::freeSpaceImpedance;
    bind_materials(m);
    bind_geometry(m);
    bind_antennas(m);
    bind_results(m);
    bind_solver(m);
}
