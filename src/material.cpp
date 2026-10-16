#include <raytracing/material.hpp>
#include <cmath>
#include <stdexcept>

namespace rt {
    Material Material::dielectric(Complex epsilon, double conductivity, double mu, std::string name) {
        return {MaterialKind::Dielectric, epsilon, mu, conductivity, std::move(name)};
    }

    Material Material::conductor(double conductivity, double mu, std::string name) {
        return {MaterialKind::Conductor, 1.0, mu, conductivity, std::move(name)};
    }

    Material Material::pec(std::string name) {
        return {MaterialKind::PEC, 1.0, 1.0, 0.0, std::move(name)};
    }

    void Material::validate() const {
        if (kind == MaterialKind::PEC) return;
        if (!(std::isfinite(epsilon.real()) && std::isfinite(epsilon.imag()) && epsilon.real() > 0 &&
              epsilon.imag() >= 0 && std::isfinite(mu) && mu > 0 && std::isfinite(conductivity) && conductivity >= 0))
            throw std::invalid_argument("Material '" + name + "' must be passive: Re(epsilon) > 0, Im(epsilon) >= 0, "
                                        "mu > 0, conductivity >= 0");
        if (kind == MaterialKind::Conductor && conductivity <= 0)
            throw std::invalid_argument("Conductor '" + name + "' needs a positive conductivity");
    }

    Complex Material::permittivity(double frequencyHz) const {
        return epsilon + Complex(0, conductivity / (2 * pi * frequencyHz * vacuumPermittivity));
    }

    Complex Material::index(double frequencyHz) const { return std::sqrt(permittivity(frequencyHz) * mu); }

    Complex Material::impedance(double frequencyHz) const {
        return freeSpaceImpedance * std::sqrt(mu / permittivity(frequencyHz));
    }

    Fresnel fresnel(const Material &incident, const Material &outgoing, double c, double frequencyHz) {
        if (outgoing.kind == MaterialKind::PEC) return {-1.0, 1.0, 0.0, 0.0};
        const Complex e1 = incident.permittivity(frequencyHz), e2 = outgoing.permittivity(frequencyHz);
        const Complex n1 = std::sqrt(e1 * incident.mu), n2 = std::sqrt(e2 * outgoing.mu);
        // Normal wavenumbers / k0; the transmitted branch decays: Im(q2) >= 0.
        const Complex q1 = n1 * c;
        Complex q2 = std::sqrt(n2 * n2 - n1 * n1 * (1 - c * c));
        if (q2.imag() < 0 || (q2.imag() == 0 && q2.real() < 0)) q2 = -q2;
        // A vanishing denominator only occurs for index matching at grazing incidence.
        const auto ratio = [](Complex a, Complex b) { return std::abs(a + b) > 0 ? (a - b) / (a + b) : Complex(0); };
        const Complex rs = ratio(outgoing.mu * q1, incident.mu * q2), rp = ratio(e2 * q1, e1 * q2);
        if (outgoing.kind == MaterialKind::Conductor) return {rs, rp, 0.0, 0.0};
        return {rs, rp, 1.0 + rs, n1 / n2 * (1.0 + rp)};
    }

    Eigen::Matrix3cd reflection_matrix(const Vec3 &incidentDirection, const Vec3 &unitNormal, Complex rs, Complex rp) {
        const Vec3 reflected = incidentDirection - 2.0 * incidentDirection.dot(unitNormal) * unitNormal;
        Vec3 s = incidentDirection.cross(unitNormal);
        // At normal incidence rp = -rs, so any perpendicular s gives the same dyadic.
        s = s.norm() > 1e-12 ? s.normalized() : unitNormal.unitOrthogonal();
        const Vec3 pIncident = s.cross(incidentDirection), pReflected = s.cross(reflected);
        return rs * (s * s.transpose()).cast<Complex>() + rp * (pReflected * pIncident.transpose()).cast<Complex>();
    }
} // namespace rt
