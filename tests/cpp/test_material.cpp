#include <raytracing/core.hpp>
#include <raytracing/solver.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
using namespace rt;
constexpr double PI = pi;

void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void near(Complex actual, Complex expected, const char *message, double tolerance = 1e-12) {
    if (std::abs(actual - expected) > tolerance) {
        std::cerr << "  actual " << actual << ", expected " << expected << '\n';
        throw std::runtime_error(message);
    }
}

template<class Action> void rejects(Action action) {
    bool rejected = false;
    try { action(); } catch (const std::invalid_argument &) { rejected = true; }
    check(rejected, "Invalid material input accepted");
}

// Textbook angle form for nonmagnetic media, written independently of the
// normal-wavenumber form used by the library (Born & Wolf; H-field p convention).
std::pair<Complex, Complex> textbook(Complex n1, Complex n2, double theta) {
    const Complex ci = std::cos(theta);
    const Complex sinT = n1 / n2 * std::sin(theta);
    // Complex - complex keeps a +0 imaginary part (double - complex would give -0,
    // selecting the growing evanescent root): principal root = decaying branch.
    const Complex ct = std::sqrt(Complex(1.0) - sinT * sinT);
    return {(n1 * ci - n2 * ct) / (n1 * ci + n2 * ct), (n2 * ci - n1 * ct) / (n2 * ci + n1 * ct)};
}
// Reflection from a half-space of material seen from the background.
struct Reflected {
    Complex s, p;
    double incidenceAngle;
};

Reflected fresnel_reflection(const Material &material, const Material &background, double c, double f) {
    const auto r = fresnel(background, material, c, f);
    return {r.rs, r.rp, std::acos(c)};
}

Material lossy(double epsilon, double lossTangent = 0, double conductivity = 0, double mu = 1) {
    return Material::dielectric(Complex(epsilon, epsilon * lossTangent), conductivity, mu);
}
} // namespace

int main() try {
    const double f = 77e9;
    const Material vacuum;

    // --- Material definitions and complex permittivity -------------------------
    const auto pec = Material::pec();
    check(pec.kind == MaterialKind::PEC && pec.name == "PEC", "Default PEC definition incorrect");
    const auto copper = Material::conductor(5.8e7);
    const double omega = 2 * PI * f;
    near(copper.permittivity(f), Complex(1, 5.8e7 / (omega * vacuumPermittivity)), "Conductivity term incorrect", 1e-6);
    near(lossy(4, 0.01).permittivity(f), Complex(4, 0.04), "Loss tangent term incorrect");
    rejects([] { lossy(4, -0.1).validate(); });
    rejects([] { Material::conductor(-1).validate(); });
    rejects([] { lossy(0).validate(); });
    rejects([] { lossy(4, 0, 0, 0).validate(); });
    near(freeSpaceImpedance, 376.730313668, "Free-space impedance constant", 1e-8);
    near(1 / std::sqrt(vacuumPermeability * vacuumPermittivity),
         speedOfLight, "Constants inconsistent", 1e-6);

    // --- PEC limit reproduces -I + 2 n n^T for every transverse field ---------
    for (const Vec3 &incident: {Vec3(1, 0, 0), Vec3(1, 1, 0).normalized(), Vec3(0.3, -0.2, 0.9).normalized(),
                                Vec3(1e-14, 0, 1).normalized(), Vec3(1, 0, 1e-7).normalized()}) {
        for (const Vec3 &normal: {Vec3(-1, 0, 0), Vec3(0, 0, 1), Vec3(1, 2, 3).normalized()}) {
            if (std::abs(incident.dot(normal)) < 1e-9) continue; // Grazing PEC is a separate case below.
            const auto coefficients = fresnel_reflection(pec, vacuum, std::abs(incident.dot(normal)), f);
            near(coefficients.s, -1, "PEC s coefficient");
            near(coefficients.p, 1, "PEC p coefficient");
            const auto R = reflection_matrix(incident, normal, coefficients.s, coefficients.p);
            const Vec3 u = incident.unitOrthogonal(), v = incident.cross(u);
            for (const Vec3C &field: {Vec3C(u.cast<Complex>()), Vec3C(v.cast<Complex>()),
                                      Vec3C((u.cast<Complex>() + Complex(0, 1) * v.cast<Complex>()) / std::sqrt(2.0))}) {
                const Vec3C n = normal.cast<Complex>();
                const Vec3C expected = -field + 2.0 * n.dot(field) * n;
                check((R * field - expected).norm() < 1e-14, "PEC dyadic differs from -I + 2nn^T");
            }
        }
    }

    // --- Lossless dielectric half-space (n = 2) -----------------------------
    const auto glass = lossy(4);
    auto normal = fresnel_reflection(glass, vacuum, 1, f);
    near(normal.s, -1.0 / 3, "Normal-incidence s coefficient");
    near(normal.p, 1.0 / 3, "Normal-incidence p coefficient (p = -s)");
    for (double degrees: {10.0, 30.0, 45.0, 60.0, 85.0}) {
        const double theta = degrees * PI / 180;
        const auto [s, p] = textbook(1.0, 2.0, theta);
        const auto r = fresnel_reflection(glass, vacuum, std::cos(theta), f);
        near(r.s, s, "Dielectric s differs from textbook Fresnel");
        near(r.p, p, "Dielectric p differs from textbook Fresnel");
        near(r.incidenceAngle, theta, "Incidence angle record");
    }
    near(fresnel_reflection(glass, vacuum, std::cos(std::atan(2.0)), f).p, 0, "Brewster angle not reproduced");
    const auto grazing = fresnel_reflection(glass, vacuum, 0, f);
    near(grazing.s, -1, "Grazing s limit");
    near(grazing.p, -1, "Grazing p limit");

    // --- Total internal reflection from a denser background ------------------
    const Material dense = lossy(4.0);
    const auto air = lossy(1);
    for (double degrees: {31.0, 45.0, 70.0, 89.0}) {
        const auto r = fresnel_reflection(air, dense, std::cos(degrees * PI / 180), f);
        near(std::abs(r.s), 1, "TIR |s| != 1");
        near(std::abs(r.p), 1, "TIR |p| != 1");
        // The principal complex root of a negative real (+0 imaginary part) is +i|.|,
        // the evanescent branch decaying into the rarer medium for exp(-i omega t).
        const auto [s, p] = textbook(2.0, 1.0, degrees * PI / 180);
        near(r.s, s, "TIR s phase incorrect");
        near(r.p, p, "TIR p phase incorrect");
    }

    // --- Impedance-matched and index-matched materials ------------------------
    near(fresnel_reflection(lossy(4, 0, 0, 4), vacuum, 1, f).s, 0,
         "Impedance-matched normal reflection");
    for (double c: {1.0, 0.5, 0.0}) {
        const auto r = fresnel_reflection(air, vacuum, c, f);
        near(r.s, 0, "Index-matched s reflection");
        near(r.p, 0, "Index-matched p reflection");
    }

    // --- Lossy media: exact normal incidence and passivity --------------------
    for (const auto &material: {copper, Material::conductor(1.4e6, 1.0),
                                lossy(5.24, 0.0, 1.38),
                                lossy(2.5, 0.3, 0.0, 2.0)}) {
        const Complex n = std::sqrt(material.permittivity(f) * material.mu);
        const Complex eta = std::sqrt(material.mu / material.permittivity(f));
        near(fresnel_reflection(material, vacuum, 1, f).s, (eta - 1.0) / (eta + 1.0),
             "Normal incidence differs from (eta2 - eta1) / (eta2 + eta1)", 1e-12);
        check(n.imag() > 0, "Lossy refractive index not decaying");
        for (int i = 0; i <= 90; ++i) {
            const auto r = fresnel_reflection(material, vacuum, std::cos(i * PI / 180), f);
            check(std::abs(r.s) <= 1 + 1e-14 && std::abs(r.p) <= 1 + 1e-14,
                  "Passive material reflects more power than incident");
            if (i < 90)
                check(std::abs(r.s) < 1 && std::abs(r.p) < 1, "Lossy reflection is lossless");
        }
    }

    // --- Good conductor: absorptivity -> 4 Rs / eta0 ---------------------------
    {
        const double rs = std::sqrt(omega * vacuumPermeability / (2 * 5.8e7));
        const double absorbed = 1 - std::norm(fresnel_reflection(copper, vacuum, 1, f).s);
        check(std::abs(absorbed / (4 * rs / freeSpaceImpedance) - 1) < 1e-3,
              "Copper absorptivity differs from surface-resistance asymptote");
        std::cout << "Copper at 77 GHz, normal incidence: 1 - |r|^2 = " << absorbed
                  << " (surface-resistance estimate " << 4 * rs / freeSpaceImpedance << ")\n";
    }

    // --- Solver: one wall reflection scales the PEC result by the coefficient --
    {
        const Rectangle wall{Vec3(1, 0, 0), Vec3::UnitY(), Vec3::UnitZ(), 5, 5};
        const auto concrete = lossy(5.24, 0.0, 1.38);
        SolverConfig config;
        config.rayCount = 4000; // Launch spacing must resolve the 0.12 m capture sphere at ~3 m.
        config.maxReflections = 1;
        const auto reflected = [](const SimulationResult &result) -> const Path & {
            for (const auto &path: result.paths)
                if (path.refinement.geometry.reflections.size() == 1) return path;
            throw std::runtime_error("Wall reflection not found");
        };
        // Tx/Rx in the xy plane: vertical isotropic antennas are s-polarized (E along z),
        // horizontal ones p-polarized (E in the plane of incidence).
        for (auto polarization: {Polarization::Vertical, Polarization::Horizontal}) {
            const Transmitter tx{Vec3::Zero(), Isotropic{polarization}};
            const Receiver rx{Vec3(0, 2, 0), Isotropic{polarization}};
            Scene pecScene, lossyScene;
            pecScene.add(wall, pec);
            lossyScene.add(wall, concrete);
            const auto a = solve(pecScene, tx, rx, config), b = solve(lossyScene, tx, rx, config);
            check(a.paths.size() == 2 && b.paths.size() == 2, "Material changed discovered geometry");
            const auto &ra = reflected(a), &rb = reflected(b);
            check(rb.field->reflections.size() == 1, "Reflection coefficients not recorded");
            near(rb.field->reflections[0].incidenceAngle, PI / 4, "Recorded incidence angle", 1e-9);
            const auto expected = fresnel_reflection(concrete, Material{}, std::cos(PI / 4), f);
            const Complex ratio = rb.coefficient / ra.coefficient;
            near(ratio, polarization == Polarization::Vertical ? -expected.s : expected.p,
                 "Solver reflection differs from the Fresnel coefficient", 1e-9);
            near(rb.field->delay, 2 * std::sqrt(2.0) / speedOfLight,
                 "Reflected delay incorrect", 1e-20);
        }
    }
    std::cout << "All material and Fresnel reflection checks passed\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
