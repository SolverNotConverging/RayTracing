#include "bindings.hpp"
#include <sstream>

void bind_materials(py::module_ &m) {
    using namespace rt;
    py::class_<Medium>(m, "Medium", "Homogeneous, lossless background medium filling the scene.")
            .def(py::init<double, double>(), py::arg("relative_permittivity") = 1.0,
                 py::arg("relative_permeability") = 1.0)
            .def_readwrite("relative_permittivity", &Medium::relativePermittivity)
            .def_readwrite("relative_permeability", &Medium::relativePermeability)
            .def("validate", &Medium::validate)
            .def("refractive_index", &Medium::refractive_index)
            .def("impedance", &Medium::impedance)
            .def("__repr__", [](const Medium &medium) {
                std::ostringstream text;
                text << "Medium(relative_permittivity=" << medium.relativePermittivity
                        << ", relative_permeability=" << medium.relativePermeability << ")";
                return text.str();
            });

    py::class_<Material>(m, "Material", R"doc(Surface material, assigned with ``Scene.add(shape, material)``.

Create one with ``Material.pec()``, ``Material.conductor(...)`` or
``Material.lossy(...)``. Each surface reflects as the planar boundary of a
half-space of its material (exact Fresnel coefficients, exp(-i omega t)):
eps = eps' (1 + i tan_delta) + i sigma / (omega eps0), evaluated at the solver
frequency. No transmitted ray is traced yet.)doc")
            .def_static("pec", &Material::pec, py::arg("name") = "PEC", "Perfect electric conductor.")
            .def_static("conductor", &Material::conductor, py::arg("name"), py::arg("conductivity"),
                        py::arg("relative_permeability") = 1.0,
                        "Good conductor from its conductivity in S/m (lattice permittivity 1).")
            .def_static("lossy", &Material::lossy, py::arg("name"), py::arg("relative_permittivity"),
                        py::arg("loss_tangent") = 0.0, py::arg("conductivity") = 0.0,
                        py::arg("relative_permeability") = 1.0,
                        "General lossy material: eps', loss tangent, conductivity (S/m) and mu_r.")
            .def_readwrite("name", &Material::name)
            .def_readwrite("perfect_conductor", &Material::perfectConductor)
            .def_readwrite("relative_permittivity", &Material::relativePermittivity)
            .def_readwrite("loss_tangent", &Material::lossTangent)
            .def_readwrite("conductivity", &Material::conductivity)
            .def_readwrite("relative_permeability", &Material::relativePermeability)
            .def("validate", &Material::validate)
            .def("complex_permittivity", &Material::complex_permittivity, py::arg("frequency_hz"),
                 "Complex relative permittivity (Im >= 0 for exp(-i omega t)).")
            .def("copy", [](const Material &material) { return material; })
            .def(py::self == py::self)
            .def("__repr__", [](const Material &material) {
                std::ostringstream text;
                if (material.perfectConductor)
                    text << "Material.pec(name='" << material.name << "')";
                else
                    text << "Material.lossy(name='" << material.name << "', relative_permittivity="
                            << material.relativePermittivity << ", loss_tangent=" << material.lossTangent
                            << ", conductivity=" << material.conductivity << ", relative_permeability="
                            << material.relativePermeability << ")";
                return text.str();
            });

    py::class_<ReflectionCoefficients>(m, "ReflectionCoefficients", R"doc(Fresnel reflection at one hit.

``s``: E perpendicular to the plane of incidence. ``p``: in-plane E, defined by
H_r = p H_i (a PEC gives s = -1, p = +1). ``incidence_angle`` is in radians.)doc")
            .def_readonly("incidence_angle", &ReflectionCoefficients::incidenceAngle)
            .def_property_readonly("incidence_angle_degrees", [](const ReflectionCoefficients &r) {
                return r.incidenceAngle * 180.0 / constants::pi;
            })
            .def_readonly("s", &ReflectionCoefficients::s)
            .def_readonly("p", &ReflectionCoefficients::p)
            .def("__repr__", [](const ReflectionCoefficients &r) {
                std::ostringstream text;
                text << "ReflectionCoefficients(incidence_angle_degrees="
                        << r.incidenceAngle * 180.0 / constants::pi << ", s=" << r.s << ", p=" << r.p << ")";
                return text.str();
            });

    m.def("fresnel_reflection", &fresnel_reflection, py::arg("material"), py::arg("background"),
          py::arg("cos_incidence"), py::arg("frequency_hz"),
          "Exact half-space Fresnel coefficients for cos(incidence angle) in [0, 1].");
    m.def("reflection_matrix", &reflection_matrix, py::arg("incident_direction"), py::arg("unit_normal"),
          py::arg("coefficients"), "3x3 complex dyadic mapping an incident transverse E to the reflected E.");
}
