#pragma once
#include "core.hpp"
#include <string>

namespace rt {
    // Dielectrics transmit; conductors and PEC only reflect.
    enum class MaterialKind { Dielectric, Conductor, PEC };

    // Passive isotropic material. The relative permittivity at frequency f is
    //   epsilon + i conductivity / (2 pi f eps0),
    // so loss appears as Im >= 0 (exp(-i omega t)). mu is real.
    struct Material {
        MaterialKind kind = MaterialKind::Dielectric;
        Complex epsilon = 1.0;
        double mu = 1.0;
        double conductivity = 0.0; // S/m
        std::string name = "vacuum";

        static Material vacuum() { return {}; }
        static Material dielectric(Complex epsilon, double conductivity = 0.0, double mu = 1.0,
                                   std::string name = "dielectric");
        static Material conductor(double conductivity, double mu = 1.0, std::string name = "conductor");
        static Material pec(std::string name = "PEC");

        void validate() const;
        Complex permittivity(double frequencyHz) const; // undefined for PEC
        Complex index(double frequencyHz) const;
        Complex impedance(double frequencyHz) const; // ohm
        bool operator==(const Material &) const = default;
    };

    // Plane-wave coefficients at a planar boundary, as functions of the incidence cosine.
    //   rs, ts: E perpendicular to the plane of incidence.
    //   rp, tp: in-plane field, defined through H (a PEC gives rs = -1, rp = +1).
    // Transmission coefficients are zero unless the outgoing material is a dielectric.
    struct Fresnel {
        Complex rs, rp, ts, tp;
    };

    Fresnel fresnel(const Material &incident, const Material &outgoing, double cosIncidence,
                    double frequencyHz);

    // 3x3 dyadic R with E_reflected = R E_incident for a transverse incident field.
    Eigen::Matrix3cd reflection_matrix(const Vec3 &incidentDirection, const Vec3 &unitNormal,
                                       Complex rs, Complex rp);
} // namespace rt
