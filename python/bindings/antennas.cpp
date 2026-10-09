#include "bindings.hpp"

void bind_antennas(py::module_ &m) {
    using namespace rt;
    py::enum_<Polarization>(m, "Polarization")
            .value("VERTICAL", Polarization::Vertical).value("HORIZONTAL", Polarization::Horizontal)
            .value("RIGHT_CIRCULAR", Polarization::RightCircular).value("LEFT_CIRCULAR", Polarization::LeftCircular);
    py::enum_<AntennaKind>(m, "AntennaKind")
            .value("ISOTROPIC", AntennaKind::Isotropic).value("SHORT_DIPOLE", AntennaKind::ShortDipole)
            .value("THIN_WIRE_DIPOLE", AntennaKind::ThinWireDipole)
            .value("RECTANGULAR_APERTURE", AntennaKind::RectangularAperture).value("IMPORTED", AntennaKind::Imported);
    py::enum_<PhasorConvention>(m, "PhasorConvention")
            .value("POSITIVE_TIME", PhasorConvention::PositiveTime).value(
                "NEGATIVE_TIME", PhasorConvention::NegativeTime);
    py::enum_<PowerReference>(m, "PowerReference")
            .value("ACCEPTED", PowerReference::Accepted).value("RADIATED", PowerReference::Radiated)
            .value("STIMULATED", PowerReference::Stimulated);
    py::class_<Medium>(m, "Medium")
            .def(py::init<double, double>(), py::arg("relative_permittivity") = 1.0,
                 py::arg("relative_permeability") = 1.0)
            .def_readwrite("relative_permittivity", &Medium::relativePermittivity)
            .def_readwrite("relative_permeability", &Medium::relativePermeability)
            .def("refractive_index", &Medium::refractive_index).def("impedance", &Medium::impedance);
    py::class_<Isotropic>(m, "Isotropic")
            .def(py::init<Polarization, Complex>(), py::arg("polarization") = Polarization::Vertical,
                 py::arg("amplitude") = Complex(1))
            .def_readwrite("polarization", &Isotropic::polarization).def_readwrite("amplitude", &Isotropic::amplitude);
    py::class_<ShortDipole>(m, "ShortDipole")
            .def(py::init<double, Complex>(), py::arg("effective_length_metres") = 1e-4,
                 py::arg("current_amperes") = Complex(1))
            .def_readwrite("effective_length_metres", &ShortDipole::effectiveLengthMetres)
            .def_readwrite("current_amperes", &ShortDipole::currentAmperes);
    py::class_<ThinWireDipole>(m, "ThinWireDipole")
            .def(py::init<double, Complex>(), py::arg("length_metres") = 0.001947,
                 py::arg("feed_current_amperes") = Complex(1))
            .def_readwrite("length_metres", &ThinWireDipole::lengthMetres)
            .def_readwrite("feed_current_amperes", &ThinWireDipole::feedCurrentAmperes);
    py::class_<RectangularAperture>(m, "RectangularAperture")
            .def(py::init<double, double, Complex, Complex, bool>(), py::arg("width_metres") = 0.01,
                 py::arg("height_metres") = 0.01, py::arg("field_x") = Complex(1), py::arg("field_y") = Complex(0),
                 py::arg("pec_backed") = true)
            .def_readwrite("width_metres", &RectangularAperture::widthMetres).def_readwrite(
                "height_metres", &RectangularAperture::heightMetres)
            .def_readwrite("field_x", &RectangularAperture::fieldX).def_readwrite(
                "field_y", &RectangularAperture::fieldY)
            .def_readwrite("pec_backed", &RectangularAperture::pecBacked);
    py::class_<FarfieldImportOptions>(m, "FarfieldImportOptions")
            .def(py::init<>()).def_readwrite("input_convention", &FarfieldImportOptions::inputConvention)
            .def_readwrite("receive_power_reference", &FarfieldImportOptions::receivePowerReference)
            .def_readwrite("input_power_watts", &FarfieldImportOptions::inputPowerWatts)
            .def_readwrite("coefficient_scale", &FarfieldImportOptions::coefficientScale)
            .def_readwrite("honor_export_axes", &FarfieldImportOptions::honorExportAxes);
    py::class_<FarfieldData, std::shared_ptr<FarfieldData> >(m, "FarfieldData")
            .def_readonly("frequencies_hz", &FarfieldData::frequenciesHz)
            .def_readonly("theta_degrees", &FarfieldData::thetaDegrees).def_readonly(
                "phi_degrees", &FarfieldData::phiDegrees)
            .def_readonly("coefficients", &FarfieldData::coefficients)
            .def_readonly("radiated_power", &FarfieldData::radiatedPower).def_readonly(
                "accepted_power", &FarfieldData::acceptedPower)
            .def_readonly("stimulated_power", &FarfieldData::stimulatedPower).def_readonly(
                "export_position", &FarfieldData::exportPosition)
            .def_readonly("export_orientation", &FarfieldData::exportOrientation).def_readonly(
                "provenance", &FarfieldData::provenance)
            .def_readonly("frequency_independent", &FarfieldData::frequencyIndependent)
            .def_readonly("input_convention", &FarfieldData::inputConvention).def_readonly(
                "coefficient_scale", &FarfieldData::coefficientScale);
    py::class_<Antenna>(m, "Antenna")
            .def(py::init<>()).def(py::init<Isotropic>(), py::arg("model"))
            .def(py::init<ShortDipole>(), py::arg("model")).def(py::init<ThinWireDipole>(), py::arg("model"))
            .def(py::init<RectangularAperture>(), py::arg("model"))
            .def_readonly("kind", &Antenna::kind).def_readwrite("orientation", &Antenna::orientation)
            .def_readwrite("isotropic", &Antenna::isotropic).def_readwrite("short_dipole", &Antenna::shortDipole)
            .def_readwrite("dipole", &Antenna::dipole).def_readwrite("aperture", &Antenna::aperture)
            .def_readwrite("receive_power_reference", &Antenna::receivePowerReference)
            .def_readwrite("receive_calibration", &Antenna::receiveCalibration)
            .def_property_readonly("pattern", [](const Antenna &a) {
                // Python exposes only read-only pattern metadata; the shared C++ data remains immutable.
                return std::const_pointer_cast<FarfieldData>(a.pattern);
            })
            .def("rotate", &Antenna::rotate, py::arg("axis"), py::arg("angle_degrees"),
                 py::return_value_policy::reference_internal)
            .def("validate", &Antenna::validate, py::arg("frequency_hz"), py::arg("medium") = Medium{})
            .def("farfield", &Antenna::farfield, py::arg("outward_direction"), py::arg("frequency_hz"),
                 py::arg("medium") = Medium{})
            .def("receive", &Antenna::receive, py::arg("look_direction"), py::arg("incident_field"),
                 py::arg("frequency_hz"), py::arg("medium") = Medium{})
            .def("reference_power", &Antenna::reference_power, py::arg("frequency_hz"), py::arg("medium") = Medium{})
            .def("copy", [](const Antenna &a) { return a; });
    py::implicitly_convertible<Isotropic, Antenna>();
    py::implicitly_convertible<ShortDipole, Antenna>();
    py::implicitly_convertible<ThinWireDipole, Antenna>();
    py::implicitly_convertible<RectangularAperture, Antenna>();
    m.def("load_farfield", &load_farfield, py::arg("path"), py::arg("options") = FarfieldImportOptions{},
          py::call_guard<py::gil_scoped_release>());
    py::class_<Transmitter>(m, "Transmitter")
            .def(py::init<Vec3, Antenna>(), py::arg("position"), py::arg("antenna") = Antenna{})
            .def_readwrite("position", &Transmitter::position).def_readwrite("antenna", &Transmitter::antenna);
    py::class_<Receiver>(m, "Receiver")
            .def(py::init<Vec3, Antenna>(), py::arg("position"), py::arg("antenna") = Antenna{})
            .def_readwrite("position", &Receiver::position).def_readwrite("antenna", &Receiver::antenna);
}
