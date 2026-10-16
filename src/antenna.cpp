#include <raytracing/antenna.hpp>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace rt {
    namespace {
        using Pattern = std::shared_ptr<const FarfieldPattern>;

        template<class... Ts>
        struct overloaded : Ts... {
            using Ts::operator()...;
        };

        bool positive(double x) { return std::isfinite(x) && x > 0; }
        bool nonzero(Complex x) { return std::isfinite(x.real()) && std::isfinite(x.imag()) && std::abs(x) > 0; }

        void require(bool condition, const char *message) {
            if (!condition) throw std::invalid_argument(message);
        }

        // Whitespace/comma separated tokens with // comments.
        struct Reader {
            std::ifstream input;
            std::istringstream row;

            explicit Reader(const std::filesystem::path &p) : input(p) {
                if (!input) throw std::runtime_error("Cannot open farfield file: " + p.string());
            }

            std::string token() {
                std::string value, line;
                while (!(row >> value)) {
                    if (!std::getline(input, line)) throw std::runtime_error("Truncated farfield file");
                    if (const auto comment = line.find("//"); comment != std::string::npos) line.resize(comment);
                    std::replace_if(line.begin(), line.end(), [](char c) { return c == ',' || c == ';'; }, ' ');
                    row.clear();
                    row.str(line);
                }
                return value;
            }

            double number() { return std::stod(token()); }
            std::size_t count() { return static_cast<std::size_t>(number()); }
            Vec3 vec() {
                const double x = number(), y = number(), z = number();
                return {x, y, z};
            }
        };

        std::string lower(std::string text) {
            std::transform(text.begin(), text.end(), text.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return text;
        }

        Vec3 theta_hat(double theta, double phi) {
            return {std::cos(theta) * std::cos(phi), std::cos(theta) * std::sin(phi), -std::sin(theta)};
        }

        Vec3 phi_hat(double phi) { return {-std::sin(phi), std::cos(phi), 0}; }

        struct Interval {
            std::size_t a, b;
            double weight;
        };

        Interval bracket(const std::vector<double> &axis, double value) {
            const double tolerance = 1e-10 * std::max({1.0, std::abs(axis.front()), std::abs(axis.back())});
            if (value < axis.front() - tolerance || value > axis.back() + tolerance)
                throw std::out_of_range("Requested direction/frequency is outside farfield coverage");
            value = std::clamp(value, axis.front(), axis.back());
            const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
            if (upper == axis.end()) return {axis.size() - 1, axis.size() - 1, 0};
            if (upper == axis.begin()) return {0, 0, 0};
            const auto b = static_cast<std::size_t>(upper - axis.begin()), a = b - 1;
            return {a, b, (value - axis[a]) / (axis[b] - axis[a])};
        }

        Interval phi_bracket(const FarfieldPattern &p, double phi) {
            const auto &axis = p.phiDegrees;
            phi = std::fmod(phi - axis.front() + 720.0, 360.0) + axis.front();
            if (phi <= axis.back()) return bracket(axis, phi);
            const double wrap = axis.front() + 360;
            if (wrap - axis.back() > 1.01 * (axis.back() - axis[axis.size() - 2]))
                throw std::out_of_range("Requested azimuth outside farfield coverage");
            return {axis.size() - 1, 0, (phi - axis.back()) / (wrap - axis.back())};
        }

        Interval frequency_bracket(const FarfieldPattern &p, double f) {
            return p.frequencyIndependent ? Interval{0, 0, 0} : bracket(p.frequenciesHz, f);
        }

        Vec3C interpolate(const FarfieldPattern &p, const Vec3 &d, double f) {
            const double theta = std::acos(std::clamp(d.z(), -1.0, 1.0)) * 180 / pi;
            const double phi = std::atan2(d.y(), d.x()) * 180 / pi;
            const auto ti = bracket(p.thetaDegrees, theta), pi_ = phi_bracket(p, phi), fi = frequency_bracket(p, f);
            const auto angular = [&](std::size_t frequency) -> Vec3C {
                const auto at = [&](std::size_t ph, std::size_t th) -> const Vec3C & {
                    return p.coefficients[(frequency * p.phiDegrees.size() + ph) * p.thetaDegrees.size() + th];
                };
                const Vec3C a = (1 - ti.weight) * at(pi_.a, ti.a) + ti.weight * at(pi_.a, ti.b);
                const Vec3C b = (1 - ti.weight) * at(pi_.b, ti.a) + ti.weight * at(pi_.b, ti.b);
                return (1 - pi_.weight) * a + pi_.weight * b;
            };
            const Vec3C value = (1 - fi.weight) * angular(fi.a) + fi.weight * angular(fi.b);
            const Vec3C dc = d.cast<Complex>();
            return value - dc.dot(value) * dc; // keep transverse after Cartesian interpolation
        }

        // Gauss-Legendre nodes and weights on [-1, 1].
        std::pair<std::vector<double>, std::vector<double>> gauss_legendre(std::size_t n) {
            std::vector<double> x(n), w(n);
            for (std::size_t i = 0; i < (n + 1) / 2; ++i) {
                double z = std::cos(pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
                double derivative = 0;
                for (int iteration = 0; iteration < 100; ++iteration) {
                    double p0 = 1, p1 = 0;
                    for (std::size_t k = 1; k <= n; ++k) {
                        const double p2 = p1;
                        p1 = p0;
                        p0 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p2) / static_cast<double>(k);
                    }
                    derivative = static_cast<double>(n) * (z * p0 - p1) / (z * z - 1.0);
                    const double step = p0 / derivative;
                    z -= step;
                    if (std::abs(step) < 1e-15) break;
                }
                x[i] = -z;
                x[n - 1 - i] = z;
                w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * derivative * derivative);
            }
            return {x, w};
        }

        double sinc(double x) { return std::abs(x) < 1e-8 ? 1 - x * x / 6 : std::sin(x) / x; }

        Vec3C isotropic_vector(const Vec3 &d, Polarization polarization) {
            const double theta = std::acos(std::clamp(d.z(), -1.0, 1.0)), phi = std::atan2(d.y(), d.x());
            const Vec3C et = theta_hat(theta, phi).cast<Complex>(), ep = phi_hat(phi).cast<Complex>();
            switch (polarization) {
                case Polarization::Vertical: return et;
                case Polarization::Horizontal: return ep;
                case Polarization::RightCircular: return (et + Complex(0, 1) * ep) / std::sqrt(2.0);
                case Polarization::LeftCircular: return (et - Complex(0, 1) * ep) / std::sqrt(2.0);
            }
            return et;
        }
    } // namespace

    Antenna &Antenna::rotate(const Vec3 &worldAxis, double angleDegrees) {
        orientation = Eigen::AngleAxisd(angleDegrees * pi / 180.0, worldAxis.normalized()).toRotationMatrix() *
                      orientation;
        return *this;
    }

    void Antenna::validate(double f, const Material &medium) const {
        require((orientation.transpose() * orientation - Eigen::Matrix3d::Identity()).norm() <= 1e-10 &&
                std::abs(orientation.determinant() - 1) <= 1e-10, "Antenna orientation must be a rotation matrix");
        const double k = wavenumber(f) * medium.index(f).real();
        std::visit(overloaded{
                       [](const Isotropic &m) { require(nonzero(m.amplitude), "Isotropic amplitude must be nonzero"); },
                       [](const ShortDipole &m) {
                           require(positive(m.effectiveLength) && nonzero(m.current),
                                   "Short dipole needs a positive length and nonzero current");
                       },
                       [k](const ThinWireDipole &m) {
                           require(positive(m.length) && nonzero(m.feedCurrent) &&
                                   std::abs(std::sin(k * m.length / 2)) >= 1e-8,
                                   "Thin-wire dipole needs a positive length, nonzero current and no feed current node");
                       },
                       [](const RectangularAperture &m) {
                           require(positive(m.width) && positive(m.height) &&
                                   std::abs(m.fieldX) + std::abs(m.fieldY) > 0,
                                   "Aperture needs positive dimensions and a nonzero field");
                       },
                       [f](const Pattern &p) {
                           require(p != nullptr, "Imported antenna has no pattern");
                           if (!p->frequencyIndependent) bracket(p->frequenciesHz, f);
                       },
                   }, model);
    }

    Vec3C Antenna::farfield(const Vec3 &world, double f, const Material &medium) const {
        const Vec3 d = orientation.transpose() * world;
        const double k = wavenumber(f) * medium.index(f).real();
        const double eta = medium.impedance(f).real();
        const Vec3 transverse = Vec3::UnitZ() - d.z() * d;
        const Vec3C local = std::visit(overloaded{
                                           [&](const Isotropic &m) -> Vec3C {
                                               return m.amplitude * isotropic_vector(d, m.polarization);
                                           },
                                           [&](const ShortDipole &m) -> Vec3C {
                                               return Complex(0, eta * k * m.effectiveLength / (4 * pi)) * m.current *
                                                      transverse.cast<Complex>();
                                           },
                                           [&](const ThinWireDipole &m) -> Vec3C {
                                               const double s = transverse.norm(), half = k * m.length / 2;
                                               if (s <= 1e-10) return Vec3C::Zero();
                                               const double shape = (std::cos(half * d.z()) - std::cos(half)) /
                                                                    (s * std::sin(half));
                                               return Complex(0, eta / (2 * pi)) * m.feedCurrent * shape *
                                                      (transverse / s).cast<Complex>();
                                           },
                                           [&](const RectangularAperture &m) -> Vec3C {
                                               if (d.z() <= 0) return Vec3C::Zero();
                                               const Vec3C e(m.fieldX, m.fieldY, 0);
                                               const Vec3C dc = d.cast<Complex>(), z = Vec3C::UnitZ();
                                               return Complex(0, -k * m.width * m.height / (2 * pi)) *
                                                      sinc(k * m.width * d.x() / 2) * sinc(k * m.height * d.y() / 2) *
                                                      dc.cross(z.cross(e));
                                           },
                                           [&](const Pattern &p) -> Vec3C { return interpolate(*p, d, f); },
                                       }, model);
        return orientation.cast<Complex>() * local;
    }

    double Antenna::reference_power(double f, const Material &medium) const {
        const double eta = medium.impedance(f).real();
        if (const auto *pattern = std::get_if<Pattern>(&model)) {
            const auto &p = **pattern;
            const auto &powers = receivePowerReference == PowerReference::Accepted ? p.acceptedPower
                                 : receivePowerReference == PowerReference::Radiated ? p.radiatedPower
                                 : p.stimulatedPower;
            const auto fi = frequency_bracket(p, f);
            const double power = (1 - fi.weight) * powers[fi.a] + fi.weight * powers[fi.b];
            require(positive(power), "Receiving with an imported pattern needs a positive power reference");
            return power;
        }
        if (const auto *m = std::get_if<Isotropic>(&model)) return 4 * pi * std::norm(m->amplitude) / (2 * eta);
        if (const auto *m = std::get_if<ShortDipole>(&model)) {
            const double k = wavenumber(f) * medium.index(f).real();
            return eta * k * k * m->effectiveLength * m->effectiveLength * std::norm(m->current) / (12 * pi);
        }
        // Integrate |rE|^2 / (2 eta) over the sphere: Gauss-Legendre in cos(theta) on
        // each hemisphere (apertures are discontinuous at the horizon) and the
        // trapezoidal rule in azimuth.
        constexpr std::size_t thetaNodes = 128, phiNodes = 256;
        const auto [nodes, weights] = gauss_legendre(thetaNodes);
        double power = 0;
        for (double hemisphere: {-1.0, 1.0})
            for (std::size_t i = 0; i < thetaNodes; ++i) {
                const double z = 0.5 * hemisphere * (nodes[i] + 1.0);
                const double rho = std::sqrt(std::max(0.0, (1.0 - z) * (1.0 + z)));
                double ring = 0;
                for (std::size_t j = 0; j < phiNodes; ++j) {
                    const double phi = 2 * pi * (static_cast<double>(j) + 0.5) / phiNodes;
                    const Vec3 local(rho * std::cos(phi), rho * std::sin(phi), z);
                    ring += farfield((orientation * local).normalized(), f, medium).squaredNorm();
                }
                power += 0.5 * weights[i] * ring * (2 * pi / phiNodes);
            }
        return power / (2 * eta);
    }

    Complex Antenna::receive(const Vec3 &look, const Vec3C &field, double f, const Material &medium) const {
        Vec3C weight;
        if (const auto *m = std::get_if<Isotropic>(&model))
            weight = orientation.cast<Complex>() * isotropic_vector(orientation.transpose() * look, m->polarization);
        else
            weight = farfield(look, f, medium) *
                     std::sqrt(4 * pi / (2 * medium.impedance(f).real() * reference_power(f, medium)));
        return receiveCalibration * (weight.transpose() * field).value();
    }

    Antenna load_farfield(const std::filesystem::path &path, const FarfieldImportOptions &options) {
        Reader r(path);
        auto p = std::make_shared<FarfieldPattern>();
        p->source = path.filename().string();
        const auto convert = [&](double theta, double phi, Complex et, Complex ep) -> Vec3C {
            const double t = theta * pi / 180, ph = phi * pi / 180;
            Vec3C value = (et * theta_hat(t, ph).cast<Complex>() + ep * phi_hat(ph).cast<Complex>()) * options.scale;
            return options.convention == PhasorConvention::PositiveTime ? value.conjugate() : value;
        };
        const auto extension = lower(path.extension().string());
        if (extension == ".ffs") {
            if (r.number() != 3.0 || r.token() != "Farfield")
                throw std::runtime_error("Only CST Farfield Source version 3.0 is supported");
            const auto nf = r.count();
            p->exportPosition = r.vec();
            const Vec3 z = r.vec(), x = r.vec();
            p->exportOrientation << x, z.cross(x), z;
            for (std::size_t i = 0; i < nf; ++i) {
                p->radiatedPower.push_back(r.number());
                p->acceptedPower.push_back(r.number());
                p->stimulatedPower.push_back(r.number());
                p->frequenciesHz.push_back(r.number());
            }
            for (std::size_t f = 0; f < nf; ++f) {
                const auto np = r.count(), nt = r.count();
                if (f == 0) {
                    p->thetaDegrees.resize(nt);
                    p->phiDegrees.resize(np);
                    p->coefficients.reserve(nf * np * nt);
                } else if (np != p->phiDegrees.size() || nt != p->thetaDegrees.size())
                    throw std::runtime_error("CST frequency blocks use different angular grids");
                for (std::size_t ph = 0; ph < np; ++ph)
                    for (std::size_t th = 0; th < nt; ++th) {
                        const double phi = r.number(), theta = r.number();
                        const double tr = r.number(), ti = r.number(), pr = r.number(), pi_ = r.number();
                        if (f == 0 && ph == 0) p->thetaDegrees[th] = theta;
                        if (f == 0 && th == 0) p->phiDegrees[ph] = phi;
                        p->coefficients.push_back(convert(theta, phi, {tr, ti}, {pr, pi_}));
                    }
            }
        } else if (extension == ".ffd") {
            const double ts = r.number(), te = r.number();
            const auto nt = r.count();
            const double ps = r.number(), pe = r.number();
            const auto np = r.count();
            for (std::size_t t = 0; t < nt; ++t) p->thetaDegrees.push_back(ts + (te - ts) * t / (nt - 1));
            for (std::size_t ph = 0; ph < np; ++ph) p->phiDegrees.push_back(ps + (pe - ps) * ph / (np - 1));
            // Either a "Frequencies n" header with Frequency blocks, or one frequency-independent block.
            auto next = r.token();
            std::size_t nf = 1;
            if (lower(next) == "frequencies") nf = r.count();
            else p->frequencyIndependent = true;
            p->coefficients.resize(nf * nt * np);
            for (std::size_t f = 0; f < nf; ++f) {
                if (p->frequencyIndependent) p->frequenciesHz.push_back(0);
                else {
                    if (lower(r.token()) != "frequency") throw std::runtime_error("Expected an HFSS Frequency block");
                    p->frequenciesHz.push_back(r.number());
                }
                for (std::size_t th = 0; th < nt; ++th)
                    for (std::size_t ph = 0; ph < np; ++ph) {
                        const double tr = p->frequencyIndependent && th == 0 && ph == 0 ? std::stod(next) : r.number();
                        const double ti = r.number(), pr = r.number(), pi_ = r.number();
                        p->coefficients[(f * np + ph) * nt + th] =
                                convert(p->thetaDegrees[th], p->phiDegrees[ph], {tr, ti}, {pr, pi_});
                    }
                p->acceptedPower.push_back(options.inputPower);
                p->radiatedPower.push_back(0);
                p->stimulatedPower.push_back(options.inputPower);
            }
        } else
            throw std::invalid_argument("Expected a CST .ffs or HFSS .ffd file");
        Antenna antenna{Pattern(p)};
        antenna.receivePowerReference = options.receivePowerReference;
        if (options.useExportAxes) antenna.orientation = p->exportOrientation;
        return antenna;
    }
} // namespace rt
