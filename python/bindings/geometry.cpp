#include "bindings.hpp"

void bind_geometry(py::module_ &m) {
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
            .def("add", &rt::Scene::add, py::arg("surface"), "Validate and copy a surface into the scene.")
            .def_property_readonly("surfaces", [](const rt::Scene &s) { return s.surfaces(); })
            .def("__len__", [](const rt::Scene &s) { return s.surfaces().size(); });
}
