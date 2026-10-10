#pragma once
#include "Surface.hpp"
#include <string>

namespace rt {
    // Homogeneous, lossless, nondispersive background medium filling the scene.
    struct Medium {
        double relativePermittivity = 1.0, relativePermeability = 1.0;

        void validate() const;

        double refractive_index() const;

        double impedance() const; // ohm
    };

    // Material assigned to a scene surface. Every surface currently acts as the
    // planar boundary of a semi-infinite half-space of its material, evaluated in
    // the local tangent plane at each hit (valid when curvature radii and the
    // surface size are large compared with the wavelength). There is no
    // transmitted ray yet: a thin dielectric sheet is modelled as a thick one.
    //
    // Phase convention exp(-i omega t): the complex relative permittivity is
    //   eps = eps' (1 + i tan_delta) + i sigma / (omega eps0),
    // so a passive material has Im(eps) >= 0. The constants are evaluated at the
    // solver frequency; dispersion beyond the conductivity term is not modelled.
    struct Material {
        std::string name = "PEC";
        bool perfectConductor = true; // The numerical parameters below are unused when true.
        double relativePermittivity = 1.0; // eps', real part of the relative permittivity
        double lossTangent = 0.0; // dielectric loss tangent, eps''/eps'
        double conductivity = 0.0; // S/m
        double relativePermeability = 1.0; // real, lossless

        static Material pec(std::string name = "PEC");

        // Good conductor described by its conductivity (lattice permittivity 1).
        static Material conductor(std::string name, double conductivity,
                                  double relativePermeability = 1.0);

        // General lossy material (dielectric, semiconductor or conductor).
        static Material lossy(std::string name, double relativePermittivity,
                              double lossTangent = 0.0, double conductivity = 0.0,
                              double relativePermeability = 1.0);

        void validate() const;

        // Complex relative permittivity at frequencyHz; throws for a PEC.
        Complex complex_permittivity(double frequencyHz) const;

        bool operator==(const Material &) const = default;
    };

    // Specular reflection coefficients at one hit.
    //   s: E perpendicular to the plane of incidence, E_r = s E_i.
    //   p: E in the plane of incidence, defined by H_r = p H_i with the tangential
    //      unit vector s_hat shared by both waves; the in-plane E unit vectors are
    //      s_hat x k_i (incident) and s_hat x k_r (reflected).
    // With this convention a PEC gives s = -1 and p = +1 at every angle, and at
    // normal incidence p = -s for any material.
    struct ReflectionCoefficients {
        double incidenceAngle = 0.0; // radians from the surface normal, in [0, pi/2]
        Complex s = -1.0, p = 1.0;
    };

    // Exact plane-wave Fresnel coefficients for a planar boundary between the
    // background medium and a half-space of material, as functions of
    // cos(incidence angle) in [0, 1].
    ReflectionCoefficients fresnel_reflection(const Material &material, const Medium &background,
                                              double cosIncidence, double frequencyHz);

    // 3x3 complex dyadic R with E_reflected = R E_incident for a transverse field.
    // incidentDirection and unitNormal must be unit vectors; the normal may face
    // either side. Near normal incidence any perpendicular s_hat gives the same R.
    Eigen::Matrix3cd reflection_matrix(const Vec3 &incidentDirection, const Vec3 &unitNormal,
                                       const ReflectionCoefficients &coefficients);
} // namespace rt
