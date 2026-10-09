#include "Antenna.hpp"
#include "EMRay.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace rt {
    namespace {
        bool positive(double x) { return std::isfinite(x) && x > 0; }

        bool finite(Complex x) {
            return std::isfinite(x.real()) && std::isfinite(x.imag());
        }

        struct Reader {
            std::ifstream input;
            std::istringstream row;

            explicit Reader(const std::filesystem::path &p) : input(p) {
                if (!input)
                    throw std::runtime_error("Cannot open farfield file: " + p.string());
            }

            std::string token() {
                std::string value, line;
                while (!(row >> value)) {
                    if (!std::getline(input, line))
                        throw std::runtime_error("Truncated farfield file");
                    if (const auto comment = line.find("//"); comment != std::string::npos)
                        line.resize(comment);
                    for (char &c: line)
                        if (c == ',' || c == ';')
                            c = ' ';
                    row.clear();
                    row.str(line);
                }
                return value;
            }

            double number() {
                const auto word = token();
                std::size_t end = 0;
                const double x = std::stod(word, &end);
                if (end != word.size() || !std::isfinite(x))
                    throw std::runtime_error("Invalid farfield number: " + word);
                return x;
            }

            std::size_t count() {
                const double n = number();
                if (n < 1 || n > 1000000 || std::floor(n) != n)
                    throw std::runtime_error("Invalid farfield sample count");
                return static_cast<std::size_t>(n);
            }

            void finish() {
                std::string extra, line;
                if (row >> extra)
                    throw std::runtime_error("Unexpected trailing farfield data");
                while (std::getline(input, line)) {
                    if (const auto c = line.find("//"); c != std::string::npos)
                        line.resize(c);
                    if (line.find_first_not_of(" \t\r\n") != std::string::npos)
                        throw std::runtime_error("Unexpected trailing farfield data");
                }
            }
        };

        Vec3 read_vec(Reader &r) {
            const double x = r.number(), y = r.number(), z = r.number();
            return {x, y, z};
        }

        Vec3C spherical(double theta, double phi, Complex et, Complex ep) {
            const double t = theta * PI / 180, p = phi * PI / 180;
            const Vec3 eTheta(std::cos(t) * std::cos(p), std::cos(t) * std::sin(p),
                              -std::sin(t));
            const Vec3 ePhi(-std::sin(p), std::cos(p), 0);
            return et * eTheta.cast<Complex>() + ep * ePhi.cast<Complex>();
        }

        struct Interval {
            std::size_t a, b;
            double weight;
        };

        Interval bracket(const std::vector<double> &axis, double value) {
            if (axis.empty())
                throw std::invalid_argument("Empty farfield axis");
            const double tolerance =
                    1e-10 * std::max({1.0, std::abs(axis.front()), std::abs(axis.back())});
            if (value < axis.front() - tolerance || value > axis.back() + tolerance)
                throw std::out_of_range(
                    "Requested direction/frequency is outside farfield coverage");
            value = std::clamp(value, axis.front(), axis.back());
            auto upper = std::upper_bound(axis.begin(), axis.end(), value);
            if (upper == axis.end())
                return {axis.size() - 1, axis.size() - 1, 0};
            if (upper == axis.begin())
                return {0, 0, 0};
            const std::size_t b = static_cast<std::size_t>(upper - axis.begin()),
                    a = b - 1;
            return {a, b, (value - axis[a]) / (axis[b] - axis[a])};
        }

        Interval phi_bracket(const FarfieldData &p, double phi) {
            const auto &axis = p.phiDegrees;
            phi = std::fmod(phi - axis.front() + 720.0, 360.0) + axis.front();
            if (phi <= axis.back())
                return bracket(axis, phi);
            const double wrap = axis.front() + 360;
            if (axis.size() < 2 ||
                wrap - axis.back() > 1.01 * (axis.back() - axis[axis.size() - 2]))
                throw std::out_of_range("Requested azimuth outside farfield coverage");
            return {axis.size() - 1, 0, (phi - axis.back()) / (wrap - axis.back())};
        }

        Vec3C interpolate(const FarfieldData &p, const Vec3 &d, double f) {
            const double theta = std::acos(std::clamp(d.z(), -1.0, 1.0)) * 180 / PI;
            const double phi = std::atan2(d.y(), d.x()) * 180 / PI;
            const auto ti = bracket(p.thetaDegrees, theta), pi = phi_bracket(p, phi);
            const auto fi =
                    p.frequencyIndependent ? Interval{0, 0, 0} : bracket(p.frequenciesHz, f);
            const auto angular = [&](std::size_t frequency) -> Vec3C {
                const auto at = [&](std::size_t ph, std::size_t th) -> const Vec3C & {
                    return p.coefficients.at(
                        (frequency * p.phiDegrees.size() + ph) * p.thetaDegrees.size() + th);
                };
                const Vec3C a =
                        (1 - ti.weight) * at(pi.a, ti.a) + ti.weight * at(pi.a, ti.b);
                const Vec3C b =
                        (1 - ti.weight) * at(pi.b, ti.a) + ti.weight * at(pi.b, ti.b);
                return (1 - pi.weight) * a + pi.weight * b;
            };
            Vec3C value = (1 - fi.weight) * angular(fi.a) + fi.weight * angular(fi.b);
            const Vec3C dc = d.cast<Complex>();
            return (value - dc.dot(value) * dc)
                    .eval(); // Preserve transversality after Cartesian interpolation.
        }

        double sinc(double x) {
            return std::abs(x) < 1e-8 ? 1 - x * x / 6 : std::sin(x) / x;
        }

        Vec3C isotropic_vector(const Vec3 &d, Polarization pol) {
            const double theta = std::acos(std::clamp(d.z(), -1.0, 1.0)),
                    phi = std::atan2(d.y(), d.x());
            const Vec3C et = Vec3(std::cos(theta) * std::cos(phi),
                                  std::cos(theta) * std::sin(phi), -std::sin(theta))
                    .cast<Complex>();
            const Vec3C ep = Vec3(-std::sin(phi), std::cos(phi), 0).cast<Complex>();
            switch (pol) {
                case Polarization::Vertical:
                    return et;
                case Polarization::Horizontal:
                    return ep;
                case Polarization::RightCircular:
                    return (et + Complex(0, 1) * ep) / std::sqrt(2.0);
                case Polarization::LeftCircular:
                    return (et - Complex(0, 1) * ep) / std::sqrt(2.0);
            }
            throw std::invalid_argument("Unknown polarization");
        }

        void ordered(const std::vector<double> &axis) {
            if (axis.empty())
                throw std::invalid_argument("Empty pattern axis");
            for (std::size_t i = 0; i < axis.size(); ++i)
                if (!std::isfinite(axis[i]) || (i && axis[i] <= axis[i - 1]))
                    throw std::invalid_argument("Pattern axes must be increasing");
        }
    } // namespace

    double Medium::refractive_index() const {
        if (!positive(relativePermittivity) || !positive(relativePermeability))
            throw std::invalid_argument(
                "Medium must have positive finite epsilon_r and mu_r");
        const double n =
                std::sqrt(relativePermittivity) * std::sqrt(relativePermeability);
        if (!positive(n))
            throw std::invalid_argument("Medium index overflow");
        return n;
    }

    double Medium::impedance() const {
        refractive_index();
        return 376.730313668 * std::sqrt(relativePermeability / relativePermittivity);
    }

    void Antenna::validate(double f, const Medium &m) const {
        if (!positive(f) || !orientation.allFinite() ||
            (orientation.transpose() * orientation - Eigen::Matrix3d::Identity())
            .norm() > 1e-10 ||
            std::abs(orientation.determinant() - 1) > 1e-10 ||
            !finite(receiveCalibration))
            throw std::invalid_argument(
                "Invalid antenna frequency, orientation or receive calibration");
        m.refractive_index();
        switch (kind) {
            case AntennaKind::Isotropic:
                if (!finite(isotropic.amplitude) || std::abs(isotropic.amplitude) == 0)
                    throw std::invalid_argument("Invalid isotropic amplitude");
                break;
            case AntennaKind::ShortDipole:
                if (!positive(shortDipole.effectiveLengthMetres) ||
                    !finite(shortDipole.currentAmperes) ||
                    std::abs(shortDipole.currentAmperes) == 0)
                    throw std::invalid_argument(
                        "Short dipole requires effective length and nonzero finite current");
                break;
            case AntennaKind::ThinWireDipole: {
                const double kl =
                        2 * PI * f * m.refractive_index() * dipole.lengthMetres / C;
                if (!positive(dipole.lengthMetres) || !finite(dipole.feedCurrentAmperes) ||
                    std::abs(dipole.feedCurrentAmperes) == 0 || !std::isfinite(kl) ||
                    std::abs(std::sin(kl / 2)) < 1e-8)
                    throw std::invalid_argument(
                        "Invalid centre-fed sinusoidal dipole (feed current node)");
                break;
            }
            case AntennaKind::RectangularAperture:
                if (!positive(aperture.widthMetres) || !positive(aperture.heightMetres) ||
                    !finite(aperture.fieldX) || !finite(aperture.fieldY) ||
                    std::abs(aperture.fieldX) + std::abs(aperture.fieldY) == 0 ||
                    !aperture.pecBacked)
                    throw std::invalid_argument("Aperture requires dimensions, complex "
                        "tangential field and PEC backing");
                break;
            case AntennaKind::Imported:
                if (!pattern)
                    throw std::invalid_argument("Missing imported farfield");
                ordered(pattern->frequenciesHz);
                ordered(pattern->thetaDegrees);
                ordered(pattern->phiDegrees);
                if (pattern->thetaDegrees.size() < 2 || pattern->phiDegrees.size() < 2 ||
                    pattern->thetaDegrees.front() < 0 ||
                    pattern->thetaDegrees.back() > 180 ||
                    pattern->coefficients.size() != pattern->frequenciesHz.size() *
                    pattern->thetaDegrees.size() *
                    pattern->phiDegrees.size())
                    throw std::invalid_argument("Invalid imported pattern grid");
                if (!pattern->frequencyIndependent)
                    bracket(pattern->frequenciesHz, f);
                break;
        }
    }

    Vec3C Antenna::farfield(const Vec3 &world, double f, const Medium &m) const {
        if (!world.allFinite() || std::abs(world.norm() - 1) > 1e-10 || !positive(f))
            throw std::invalid_argument(
                "Antenna evaluation requires a unit direction and positive frequency");
        const Vec3 d = orientation.transpose() * world;
        const double k = 2 * PI * f * m.refractive_index() / C;
        Vec3C value = Vec3C::Zero();
        const Vec3 transverse = Vec3::UnitZ() - d.z() * d;
        switch (kind) {
            case AntennaKind::Isotropic:
                value = isotropic.amplitude * isotropic_vector(d, isotropic.polarization);
                break;
            case AntennaKind::ShortDipole:
                value = Complex(0, m.impedance() * k * shortDipole.effectiveLengthMetres /
                                   (4 * PI)) *
                        shortDipole.currentAmperes * transverse.cast<Complex>();
                break;
            case AntennaKind::ThinWireDipole: {
                const double s = transverse.norm(), half = k * dipole.lengthMetres / 2;
                if (s > 1e-10) {
                    const double shape =
                            (std::cos(half * d.z()) - std::cos(half)) / (s * std::sin(half));
                    value = Complex(0, m.impedance() / (2 * PI)) * dipole.feedCurrentAmperes *
                            shape * (transverse / s).cast<Complex>();
                }
                break;
            }
            case AntennaKind::RectangularAperture:
                if (d.z() > 0) {
                    const Vec3C e(aperture.fieldX, aperture.fieldY, 0);
                    const Vec3C dc = d.cast<Complex>(), z = Vec3::UnitZ().cast<Complex>();
                    value = Complex(0, -k * aperture.widthMetres * aperture.heightMetres /
                                       (2 * PI)) *
                            sinc(k * aperture.widthMetres * d.x() / 2) *
                            sinc(k * aperture.heightMetres * d.y() / 2) *
                            dc.cross(z.cross(e));
                }
                break;
            case AntennaKind::Imported:
                if (!pattern)
                    throw std::invalid_argument("Missing farfield samples");
                value = interpolate(*pattern, d, f);
                break;
        }
        const Vec3C result = orientation.cast<Complex>() * value;
        if (!result.allFinite())
            throw std::invalid_argument("Nonfinite antenna field");
        return result;
    }

    double Antenna::reference_power(double f, const Medium &m) const {
        if (kind == AntennaKind::Imported) {
            if (!pattern)
                throw std::invalid_argument("Missing imported farfield");
            const auto &powers = receivePowerReference == PowerReference::Accepted
                                     ? pattern->acceptedPower
                                     : receivePowerReference == PowerReference::Radiated
                                           ? pattern->radiatedPower
                                           : pattern->stimulatedPower;
            if (powers.size() != pattern->frequenciesHz.size())
                throw std::invalid_argument("Missing receive input-power metadata");
            const auto fi = pattern->frequencyIndependent
                                ? Interval{0, 0, 0}
                                : bracket(pattern->frequenciesHz, f);
            const double power =
                    (1 - fi.weight) * powers[fi.a] + fi.weight * powers[fi.b];
            if (!positive(power))
                throw std::invalid_argument(
                    "Receive farfield requires an explicit positive power reference");
            return power;
        }
        if (kind == AntennaKind::Isotropic)
            return 4 * PI * std::norm(isotropic.amplitude) / (2 * m.impedance());
        if (kind == AntennaKind::ShortDipole) {
            const double k = 2 * PI * f * m.refractive_index() / C;
            return m.impedance() * k * k * shortDipole.effectiveLengthMetres *
                   shortDipole.effectiveLengthMetres *
                   std::norm(shortDipole.currentAmperes) / (12 * PI);
        }
        double power = 0;
        constexpr std::size_t count = 4096;
        for (const auto &d: launch_directions(count))
            power += farfield(d, f, m).squaredNorm();
        power *= 4 * PI / (count * 2 * m.impedance());
        if (!positive(power))
            throw std::invalid_argument("Antenna has zero radiated power");
        return power;
    }

    Complex Antenna::receive(const Vec3 &look, const Vec3C &field, double f,
                             const Medium &m) const {
        if (!field.allFinite())
            throw std::invalid_argument("Nonfinite incident field");
        Vec3C weight;
        if (kind == AntennaKind::Isotropic)
            weight = orientation.cast<Complex>() *
                     isotropic_vector(orientation.transpose() * look,
                                      isotropic.polarization);
        else
            weight = farfield(look, f, m) *
                     std::sqrt(4 * PI / (2 * m.impedance() * reference_power(f, m)));
        return receiveCalibration * (weight.transpose() * field)(0, 0);
    }

    Antenna load_farfield(const std::filesystem::path &path,
                          const FarfieldImportOptions &options) {
        if (!positive(options.coefficientScale) ||
            !std::isfinite(options.inputPowerWatts) || options.inputPowerWatts < 0)
            throw std::invalid_argument("Invalid farfield scale or power reference");
        Reader r(path);
        auto p = std::make_shared<FarfieldData>();
        p->provenance = path.filename().string();
        p->inputConvention = options.inputConvention;
        p->coefficientScale = options.coefficientScale;
        const auto convert = [&](double theta, double phi, Complex et,
                                 Complex ep) -> Vec3C {
            Vec3C value = spherical(theta, phi, et, ep) * options.coefficientScale;
            if (options.inputConvention == PhasorConvention::PositiveTime)
                value = value.conjugate().eval();
            return value;
        };
        auto ext = path.extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        if (ext == ".ffs") {
            if (std::abs(r.number() - 3.0) > 1e-12 || r.token() != "Farfield")
                throw std::runtime_error(
                    "Only CST Farfield Source version 3.0 is supported");
            const auto nf = r.count();
            p->exportPosition = read_vec(r);
            const Vec3 z = read_vec(r), x = read_vec(r);
            p->exportOrientation.col(0) = x;
            p->exportOrientation.col(1) = z.cross(x);
            p->exportOrientation.col(2) = z;
            for (std::size_t i = 0; i < nf; ++i) {
                p->radiatedPower.push_back(r.number());
                p->acceptedPower.push_back(r.number());
                p->stimulatedPower.push_back(r.number());
                p->frequenciesHz.push_back(r.number());
            }
            for (std::size_t f = 0; f < nf; ++f) {
                const auto np = r.count(), nt = r.count();
                if (np < 2 || nt < 2 || np * nt > 10000000 || nf * np * nt > 20000000)
                    throw std::runtime_error("Unsupported farfield grid size");
                if (f == 0) {
                    p->thetaDegrees.resize(nt);
                    p->phiDegrees.resize(np);
                    p->coefficients.reserve(nf * np * nt);
                } else if (np != p->phiDegrees.size() || nt != p->thetaDegrees.size())
                    throw std::runtime_error("Frequency grids differ");
                for (std::size_t ph = 0; ph < np; ++ph)
                    for (std::size_t th = 0; th < nt; ++th) {
                        const double phi = r.number(), theta = r.number();
                        const double tr = r.number(), ti = r.number(), pr = r.number(),
                                pi = r.number();
                        if (f == 0 && ph == 0)
                            p->thetaDegrees[th] = theta;
                        if (f == 0 && th == 0)
                            p->phiDegrees[ph] = phi;
                        if (std::abs(p->thetaDegrees[th] - theta) > 1e-8 ||
                            std::abs(p->phiDegrees[ph] - phi) > 1e-8)
                            throw std::runtime_error("CST angular grid is inconsistent");
                        p->coefficients.push_back(convert(theta, phi, {tr, ti}, {pr, pi}));
                    }
            }
        } else if (ext == ".ffd") {
            const double ts = r.number(), te = r.number();
            const auto nt = r.count();
            const double ps = r.number(), pe = r.number();
            const auto np = r.count();
            if (nt < 2 || np < 2 || nt * np > 10000000)
                throw std::runtime_error("Invalid HFSS angular grid");
            for (std::size_t t = 0; t < nt; ++t)
                p->thetaDegrees.push_back(ts + (te - ts) * t / (nt - 1));
            for (std::size_t ph = 0; ph < np; ++ph)
                p->phiDegrees.push_back(ps + (pe - ps) * ph / (np - 1));
            auto next = r.token();
            auto lower = next;
            std::transform(
                lower.begin(), lower.end(), lower.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            std::size_t nf = 1;
            if (lower == "frequencies")
                nf = r.count();
            else
                p->frequencyIndependent = true;
            if (nf * nt * np > 20000000)
                throw std::runtime_error("Farfield grid too large");
            p->coefficients.resize(nf * nt * np);
            for (std::size_t f = 0; f < nf; ++f) {
                if (!p->frequencyIndependent) {
                    auto label = r.token();
                    std::transform(
                        label.begin(), label.end(), label.begin(),
                        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                    if (label != "frequency")
                        throw std::runtime_error("Expected HFSS Frequency block");
                    p->frequenciesHz.push_back(r.number());
                } else
                    p->frequenciesHz.push_back(0);
                for (std::size_t th = 0; th < nt; ++th)
                    for (std::size_t ph = 0; ph < np; ++ph) {
                        double tr;
                        if (p->frequencyIndependent && th == 0 && ph == 0) {
                            std::size_t used = 0;
                            tr = std::stod(next, &used);
                            if (used != next.size() || !std::isfinite(tr))
                                throw std::runtime_error("Invalid first HFSS field sample");
                        } else
                            tr = r.number();
                        const double ti = r.number(), pr = r.number(), pi = r.number();
                        p->coefficients[(f * np + ph) * nt + th] = convert(
                            p->thetaDegrees[th], p->phiDegrees[ph], {tr, ti}, {pr, pi});
                    }
                p->acceptedPower.push_back(options.inputPowerWatts);
                p->radiatedPower.push_back(0);
                p->stimulatedPower.push_back(options.inputPowerWatts);
            }
        } else
            throw std::invalid_argument("Expected CST .ffs or HFSS .ffd file");
        r.finish();
        Antenna antenna;
        antenna.kind = AntennaKind::Imported;
        antenna.pattern = p;
        antenna.receivePowerReference = options.receivePowerReference;
        if (options.honorExportAxes)
            antenna.orientation = p->exportOrientation;
        antenna.validate(p->frequencyIndependent ? 1.0 : p->frequenciesHz.front());
        for (const auto &value: p->coefficients)
            if (!value.allFinite())
                throw std::runtime_error("Invalid farfield coefficient");
        return antenna;
    }
} // namespace rt
