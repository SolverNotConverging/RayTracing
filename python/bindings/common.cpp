// Constants, geometry, materials and antennas shared by both solvers.
#include "bindings.hpp"

using namespace rt;

namespace {
    using Pattern = std::shared_ptr<const FarfieldPattern>;
} // namespace

void bind_common(py::module_ &m) {
    m.attr("SPEED_OF_LIGHT") = speedOfLight;
    m.attr("VACUUM_PERMEABILITY") = vacuumPermeability;
    m.attr("VACUUM_PERMITTIVITY") = vacuumPermittivity;
    m.attr("FREE_SPACE_IMPEDANCE") = freeSpaceImpedance;

    // ---- Geometry
    py::class_<Rectangle>(m, "Rectangle", "Flat rectangle; u and v are perpendicular unit axes, normal u x v.")
            .def(py::init<Vec3, Vec3, Vec3, double, double>(), "center"_a, "u"_a, "v"_a, "half_width"_a,
                 "half_height"_a)
            .def_readwrite("center", &Rectangle::center).def_readwrite("u", &Rectangle::u)
            .def_readwrite("v", &Rectangle::v).def_readwrite("half_width", &Rectangle::halfWidth)
            .def_readwrite("half_height", &Rectangle::halfHeight);
    py::class_<Disk>(m, "Disk")
            .def(py::init<Vec3, Vec3, double>(), "center"_a, "normal"_a, "radius"_a)
            .def_readwrite("center", &Disk::center).def_readwrite("normal", &Disk::normal)
            .def_readwrite("radius", &Disk::radius);
    py::class_<Sphere>(m, "Sphere")
            .def(py::init<Vec3, double>(), "center"_a, "radius"_a)
            .def_readwrite("center", &Sphere::center).def_readwrite("radius", &Sphere::radius);
    py::class_<Cylinder>(m, "Cylinder", "Finite cylinder; center is the axis midpoint.")
            .def(py::init<Vec3, Vec3, double, double, bool>(), "center"_a, "axis"_a, "radius"_a, "half_length"_a,
                 "capped"_a = true)
            .def_readwrite("center", &Cylinder::center).def_readwrite("axis", &Cylinder::axis)
            .def_readwrite("radius", &Cylinder::radius).def_readwrite("half_length", &Cylinder::halfLength)
            .def_readwrite("capped", &Cylinder::capped);
    py::class_<Triangle>(m, "Triangle")
            .def(py::init<Vec3, Vec3, Vec3>(), "a"_a, "b"_a, "c"_a)
            .def_readwrite("a", &Triangle::a).def_readwrite("b", &Triangle::b).def_readwrite("c", &Triangle::c);
    py::class_<Box>(m, "Box", "Closed box; axes columns form a right-handed frame.")
            .def(py::init<Vec3, Vec3, Eigen::Matrix3d>(), "center"_a, "half_lengths"_a,
                 "axes"_a = Eigen::Matrix3d::Identity())
            .def_readwrite("center", &Box::center).def_readwrite("half_lengths", &Box::halfLengths)
            .def_readwrite("axes", &Box::axes);
    m.def("faces", &faces, "shape"_a, "Surfaces bounding a shape (six rectangles for a Box).");

    // ---- Materials
    py::enum_<MaterialKind>(m, "MaterialKind")
            .value("DIELECTRIC", MaterialKind::Dielectric).value("CONDUCTOR", MaterialKind::Conductor)
            .value("PEC", MaterialKind::PEC);
    py::class_<Material>(m, "Material", "Passive isotropic material; permittivity is epsilon + i sigma/(omega eps0).")
            .def(py::init<>())
            .def_static("vacuum", &Material::vacuum)
            .def_static("dielectric", &Material::dielectric, "epsilon"_a, "conductivity"_a = 0.0, "mu"_a = 1.0,
                        "name"_a = "dielectric")
            .def_static("conductor", &Material::conductor, "conductivity"_a, "mu"_a = 1.0, "name"_a = "conductor")
            .def_static("pec", &Material::pec, "name"_a = "PEC")
            .def_readonly("kind", &Material::kind).def_readonly("epsilon", &Material::epsilon)
            .def_readonly("mu", &Material::mu).def_readonly("conductivity", &Material::conductivity)
            .def_readonly("name", &Material::name)
            .def("permittivity", &Material::permittivity, "frequency_hz"_a)
            .def("index", &Material::index, "frequency_hz"_a)
            .def("impedance", &Material::impedance, "frequency_hz"_a)
            .def(py::self == py::self)
            .def("__repr__", [](const Material &x) {
                return "Material(" + x.name + ", " + py::repr(py::cast(x.kind)).cast<std::string>() + ", epsilon=" +
                       py::repr(py::cast(x.epsilon)).cast<std::string>() + ", mu=" + std::to_string(x.mu) +
                       ", conductivity=" + std::to_string(x.conductivity) + ")";
            });
    py::class_<Fresnel>(m, "Fresnel")
            .def_readonly("rs", &Fresnel::rs).def_readonly("rp", &Fresnel::rp)
            .def_readonly("ts", &Fresnel::ts).def_readonly("tp", &Fresnel::tp);
    m.def("fresnel", &fresnel, "incident"_a, "outgoing"_a, "cos_incidence"_a, "frequency_hz"_a);

    // ---- Antennas
    py::enum_<Polarization>(m, "Polarization")
            .value("VERTICAL", Polarization::Vertical).value("HORIZONTAL", Polarization::Horizontal)
            .value("RIGHT_CIRCULAR", Polarization::RightCircular).value("LEFT_CIRCULAR", Polarization::LeftCircular);
    py::enum_<PhasorConvention>(m, "PhasorConvention")
            .value("POSITIVE_TIME", PhasorConvention::PositiveTime)
            .value("NEGATIVE_TIME", PhasorConvention::NegativeTime);
    py::enum_<PowerReference>(m, "PowerReference")
            .value("ACCEPTED", PowerReference::Accepted).value("RADIATED", PowerReference::Radiated)
            .value("STIMULATED", PowerReference::Stimulated);
    py::class_<Isotropic>(m, "Isotropic")
            .def(py::init<Polarization, Complex>(), "polarization"_a = Polarization::Vertical, "amplitude"_a = 1.0)
            .def_readwrite("polarization", &Isotropic::polarization)
            .def_readwrite("amplitude", &Isotropic::amplitude);
    py::class_<ShortDipole>(m, "ShortDipole", "Uniform-current element along local z.")
            .def(py::init<double, Complex>(), "effective_length"_a = 1e-4, "current"_a = 1.0)
            .def_readwrite("effective_length", &ShortDipole::effectiveLength)
            .def_readwrite("current", &ShortDipole::current);
    py::class_<ThinWireDipole>(m, "ThinWireDipole", "Centre-fed sinusoidal-current wire along local z.")
            .def(py::init<double, Complex>(), "length"_a = 0.001947, "feed_current"_a = 1.0)
            .def_readwrite("length", &ThinWireDipole::length)
            .def_readwrite("feed_current", &ThinWireDipole::feedCurrent);
    py::class_<RectangularAperture>(m, "RectangularAperture", "PEC-backed aperture radiating into local +z.")
            .def(py::init<double, double, Complex, Complex>(), "width"_a = 0.01, "height"_a = 0.01,
                 "field_x"_a = 1.0, "field_y"_a = 0.0)
            .def_readwrite("width", &RectangularAperture::width)
            .def_readwrite("height", &RectangularAperture::height)
            .def_readwrite("field_x", &RectangularAperture::fieldX)
            .def_readwrite("field_y", &RectangularAperture::fieldY);
    py::class_<FarfieldPattern, std::shared_ptr<FarfieldPattern>>(m, "FarfieldPattern")
            .def_readonly("frequencies_hz", &FarfieldPattern::frequenciesHz)
            .def_readonly("theta_degrees", &FarfieldPattern::thetaDegrees)
            .def_readonly("phi_degrees", &FarfieldPattern::phiDegrees)
            .def_readonly("coefficients", &FarfieldPattern::coefficients)
            .def_readonly("radiated_power", &FarfieldPattern::radiatedPower)
            .def_readonly("accepted_power", &FarfieldPattern::acceptedPower)
            .def_readonly("stimulated_power", &FarfieldPattern::stimulatedPower)
            .def_readonly("export_position", &FarfieldPattern::exportPosition)
            .def_readonly("export_orientation", &FarfieldPattern::exportOrientation)
            .def_readonly("frequency_independent", &FarfieldPattern::frequencyIndependent)
            .def_readonly("source", &FarfieldPattern::source);
    py::class_<FarfieldImportOptions>(m, "FarfieldImportOptions")
            .def(py::init([](PhasorConvention convention, PowerReference reference, double inputPower, double scale,
                             bool useExportAxes) {
                     return FarfieldImportOptions{convention, reference, inputPower, scale, useExportAxes};
                 }), "convention"_a = PhasorConvention::PositiveTime, "receive_power_reference"_a = PowerReference::Accepted,
                 "input_power"_a = 0.0, "scale"_a = 1.0, "use_export_axes"_a = true)
            .def_readwrite("convention", &FarfieldImportOptions::convention)
            .def_readwrite("receive_power_reference", &FarfieldImportOptions::receivePowerReference)
            .def_readwrite("input_power", &FarfieldImportOptions::inputPower)
            .def_readwrite("scale", &FarfieldImportOptions::scale)
            .def_readwrite("use_export_axes", &FarfieldImportOptions::useExportAxes);

    py::class_<Antenna>(m, "Antenna")
            .def(py::init<>())
            .def(py::init<Isotropic>(), "model"_a).def(py::init<ShortDipole>(), "model"_a)
            .def(py::init<ThinWireDipole>(), "model"_a).def(py::init<RectangularAperture>(), "model"_a)
            .def_property_readonly("model", [](const Antenna &a) -> py::object {
                return std::visit([](const auto &model) -> py::object {
                    if constexpr (std::is_same_v<std::decay_t<decltype(model)>, Pattern>)
                        return py::cast(std::const_pointer_cast<FarfieldPattern>(model));
                    else return py::cast(model);
                }, a.model);
            })
            .def_readwrite("orientation", &Antenna::orientation)
            .def_readwrite("receive_power_reference", &Antenna::receivePowerReference)
            .def_readwrite("receive_calibration", &Antenna::receiveCalibration)
            .def("rotate", &Antenna::rotate, "axis"_a, "angle_degrees"_a, py::return_value_policy::reference_internal)
            .def("validate", &Antenna::validate, "frequency_hz"_a, "medium"_a = Material{})
            .def("farfield", &Antenna::farfield, "direction"_a, "frequency_hz"_a, "medium"_a = Material{})
            .def("receive", &Antenna::receive, "look_direction"_a, "incident_field"_a, "frequency_hz"_a,
                 "medium"_a = Material{})
            .def("reference_power", &Antenna::reference_power, "frequency_hz"_a, "medium"_a = Material{});
    py::implicitly_convertible<Isotropic, Antenna>();
    py::implicitly_convertible<ShortDipole, Antenna>();
    py::implicitly_convertible<ThinWireDipole, Antenna>();
    py::implicitly_convertible<RectangularAperture, Antenna>();
    m.def("load_farfield", &load_farfield, "path"_a, "options"_a = FarfieldImportOptions{});
}
