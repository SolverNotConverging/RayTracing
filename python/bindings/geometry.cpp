#include "bindings.hpp"

void bind_geometry(py::module_ &m) {
    m.def("_validate_surface", &validate_surface);
    m.def("_intersect", [](const Surface &shape, const Vec3 &origin, const Vec3 &direction,
                           double tMin, double tMax) -> py::object {
        const auto hit = std::visit([&](const auto &s) { return intersect(origin, direction, s, tMin, tMax); }, shape);
        if (!hit) return py::none();
        return py::make_tuple(hit->distance_, hit->normal_);
    });
    py::class_<Rectangle>(m, "Rectangle")
            .def(py::init<Vec3, Vec3, Vec3, double, double>(), py::arg("center"), py::arg("u"),
                 py::arg("v"), py::arg("half_width"), py::arg("half_height"))
            .def_readwrite("center", &Rectangle::center_)
            .def_readwrite("u", &Rectangle::u_).def_readwrite("v", &Rectangle::v_)
            .def_readwrite("half_width", &Rectangle::halfWidth_)
            .def_readwrite("half_height", &Rectangle::halfHeight_);
    py::class_<Disk>(m, "Disk")
            .def(py::init<Vec3, Vec3, double>(), py::arg("center"), py::arg("normal"), py::arg("radius"))
            .def_readwrite("center", &Disk::center_).def_readwrite("normal", &Disk::normal_)
            .def_readwrite("radius", &Disk::radius_);
    py::class_<Sphere>(m, "Sphere")
            .def(py::init<Vec3, double>(), py::arg("center"), py::arg("radius"))
            .def_readwrite("center", &Sphere::center_).def_readwrite("radius", &Sphere::radius_);
    py::class_<Cylinder>(m, "Cylinder", "Finite cylinder: center is the axis midpoint; lengths are in metres.")
            .def(py::init<Vec3, Vec3, double, double, bool>(), py::arg("center"), py::arg("axis"),
                 py::arg("radius"), py::arg("half_length"), py::arg("capped") = true)
            .def_readwrite("center", &Cylinder::center_).def_readwrite("axis", &Cylinder::axis_)
            .def_readwrite("radius", &Cylinder::radius_).def_readwrite("half_length", &Cylinder::halfLength_)
            .def_readwrite("capped", &Cylinder::capped_);
    py::class_<Triangle>(m, "Triangle")
            .def(py::init<Vec3, Vec3, Vec3>(), py::arg("a"), py::arg("b"), py::arg("c"))
            .def_readwrite("a", &Triangle::a_).def_readwrite("b", &Triangle::b_).def_readwrite("c", &Triangle::c_);
    py::class_<rt::Scene>(m, "Scene")
            .def(py::init<>())
            .def("add", &rt::Scene::add, py::arg("surface"), py::arg("material"),
                 "Validate and copy a surface and its material into the scene; return the surface index.")
            .def_property_readonly("surfaces", [](const rt::Scene &s) { return s.surfaces(); })
            .def_property_readonly("materials", [](const rt::Scene &s) { return s.materials(); },
                                   "Material of each surface, indexed like surfaces.")
            .def("__len__", [](const rt::Scene &s) { return s.surfaces().size(); });
}
