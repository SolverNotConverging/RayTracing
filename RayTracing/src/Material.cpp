#include "Material.hpp"
#include "Constants.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rt {
    namespace {
        bool positive(double x) { return std::isfinite(x) && x > 0.0; }

        bool nonnegative(double x) { return std::isfinite(x) && x >= 0.0; }

        void require_unit(const Vec3 &v, const char *message) {
            if (!v.allFinite() || std::abs(v.norm() - 1.0) > 1e-10)
                throw std::invalid_argument(message);
        }

        // (a q1 - b q2) / (a q1 + b q2). A vanishing denominator only occurs for
        // an index-matched lossless half-space at grazing incidence, where the
        // limit of the reflection coefficient is zero.
        Complex ratio(Complex a, Complex q1, Complex b, Complex q2) {
            const Complex denominator = a * q1 + b * q2;
            if (denominator == Complex(0.0))
                return 0.0;
            return (a * q1 - b * q2) / denominator;
        }
    } // namespace

    void Medium::validate() const {
        if (!positive(relativePermittivity) || !positive(relativePermeability))
            throw std::invalid_argument("Medium must have positive finite epsilon_r and mu_r");
    }

    double Medium::refractive_index() const {
        validate();
        const double n = std::sqrt(relativePermittivity) * std::sqrt(relativePermeability);
        if (!positive(n))
            throw std::invalid_argument("Medium index overflow");
        return n;
    }

    double Medium::impedance() const {
        validate();
        return constants::freeSpaceImpedance * std::sqrt(relativePermeability / relativePermittivity);
    }

    Material Material::pec(std::string name) {
        Material material;
        material.name = std::move(name);
        return material;
    }

    Material Material::conductor(std::string name, double conductivity, double relativePermeability) {
        return lossy(std::move(name), 1.0, 0.0, conductivity, relativePermeability);
    }

    Material Material::lossy(std::string name, double relativePermittivity, double lossTangent,
                             double conductivity, double relativePermeability) {
        Material material;
        material.name = std::move(name);
        material.perfectConductor = false;
        material.relativePermittivity = relativePermittivity;
        material.lossTangent = lossTangent;
        material.conductivity = conductivity;
        material.relativePermeability = relativePermeability;
        material.validate();
        return material;
    }

    void Material::validate() const {
        if (perfectConductor)
            return;
        if (!positive(relativePermittivity) || !nonnegative(lossTangent) ||
            !nonnegative(conductivity) || !positive(relativePermeability))
            throw std::invalid_argument(
                "Material '" + name +
                "' requires positive eps' and mu_r, and nonnegative loss tangent and conductivity");
    }

    Complex Material::complex_permittivity(double frequencyHz) const {
        if (perfectConductor)
            throw std::invalid_argument("A perfect conductor has no finite permittivity");
        validate();
        if (!positive(frequencyHz))
            throw std::invalid_argument("Material frequency must be positive and finite");
        const double omega = 2.0 * constants::pi * frequencyHz;
        const Complex epsilon(relativePermittivity,
                              relativePermittivity * lossTangent +
                              conductivity / (omega * constants::vacuumPermittivity));
        if (!std::isfinite(epsilon.real()) || !std::isfinite(epsilon.imag()))
            throw std::invalid_argument("Material permittivity overflow");
        return epsilon;
    }

    ReflectionCoefficients fresnel_reflection(const Material &material, const Medium &background,
                                              double cosIncidence, double frequencyHz) {
        material.validate();
        background.validate();
        if (!std::isfinite(cosIncidence) || cosIncidence < -1e-12 || cosIncidence > 1.0 + 1e-12)
            throw std::invalid_argument("cos(incidence angle) must lie in [0, 1]");
        if (!positive(frequencyHz))
            throw std::invalid_argument("Reflection frequency must be positive and finite");
        const double c = std::clamp(cosIncidence, 0.0, 1.0);
        ReflectionCoefficients result;
        result.incidenceAngle = std::acos(c);
        if (material.perfectConductor)
            return result; // s = -1, p = +1

        const double e1 = background.relativePermittivity, m1 = background.relativePermeability;
        const Complex e2 = material.complex_permittivity(frequencyHz);
        const double m2 = material.relativePermeability;
        const double sin2 = (1.0 - c) * (1.0 + c);
        // Normalized normal wavenumbers k_z / k0 on each side. The transmitted
        // branch decays away from the boundary, Im(q2) >= 0 for exp(-i omega t);
        // Im(e2) >= 0 makes the principal square root that branch, provided a
        // zero imaginary part is +0 (never -0).
        const double q1 = std::sqrt(e1 * m1) * c;
        Complex argument = e2 * m2 - e1 * m1 * sin2;
        if (!(argument.imag() > 0.0))
            argument = {argument.real(), 0.0};
        const Complex q2 = std::sqrt(argument);
        result.s = ratio(m2, q1, m1, q2);
        result.p = ratio(e2, q1, e1, q2);
        if (!std::isfinite(result.s.real()) || !std::isfinite(result.s.imag()) ||
            !std::isfinite(result.p.real()) || !std::isfinite(result.p.imag()))
            throw std::invalid_argument("Nonfinite Fresnel coefficient");
        return result;
    }

    Eigen::Matrix3cd reflection_matrix(const Vec3 &incidentDirection, const Vec3 &unitNormal,
                                       const ReflectionCoefficients &coefficients) {
        require_unit(incidentDirection, "Reflection requires a unit incident direction");
        require_unit(unitNormal, "Reflection requires a unit surface normal");
        const Vec3 reflected = incidentDirection - 2.0 * incidentDirection.dot(unitNormal) * unitNormal;
        Vec3 s = incidentDirection.cross(unitNormal);
        // At normal incidence the plane of incidence is undefined; R reduces to
        // s * (transverse identity) because p = -s, so any perpendicular works.
        s = s.norm() > 1e-12 ? s.normalized() : unitNormal.unitOrthogonal();
        const Vec3 pIncident = s.cross(incidentDirection);
        const Vec3 pReflected = s.cross(reflected);
        return coefficients.s * (s * s.transpose()).cast<Complex>() +
               coefficients.p * (pReflected * pIncident.transpose()).cast<Complex>();
    }
} // namespace rt
