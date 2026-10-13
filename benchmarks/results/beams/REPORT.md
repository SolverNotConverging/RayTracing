# Beam transport analytical results

Solver 0.5.0; all checks passed: True.

API: `raytracing.transport.solve`; finite Gaussian beams and plane-wave interface limits.

| Check | Maximum error | Limit | Pass |
| --- | ---: | ---: | --- |
| Gaussian complex field across waist | 5.882e-14 | 2e-11 | True |
| Gaussian receiver quadrature | 3.179e-05 | 4e-05 | True |
| Gaussian quadrature error reduction | 0.000222 | 0.001 | True |
| Slab coherent amplitude eps=(4+0j) | 4.334e-15 | 2e-11 | True |
| Slab convergence eps=(4+0j) | 3.901e-14 | 1e-08 | True |
| Slab coherent amplitude eps=(4+0.02j) | 4.314e-15 | 2e-11 | True |
| Slab convergence eps=(4+0.02j) | 3.99e-14 | 1e-08 | True |
| Slab branch energy accounting | 3.268e-16 | 2e-12 | True |
| Finite Gaussian slab coherent transmission | 4.055e-15 | 2e-11 | True |
| Pruning convergence | 4.334e-15 | 2e-11 | True |
| Fresnel angle sweep | 3.331e-16 | 2e-12 | True |
| Snell refraction | 5.551e-17 | 2e-12 | True |
| Lossless R+T | 4.441e-16 | 2e-12 | True |
| Concave mirror complex field and Gouy phase | 1.142e-13 | 2e-11 | True |
| Concave mirror beam width | 5.551e-17 | 2e-12 | True |
| Oblique spherical mirror astigmatism | 4.441e-16 | 2e-11 | True |
| Repeated spherical/cylindrical Gaussian focus complex field | 1.556e-15 | 2e-10 | True |
| Repeated spherical/cylindrical Gaussian focus widths | 8.882e-16 | 2e-11 | True |
| Gaussian curved refraction complex field | 8.528e-14 | 2e-11 | True |
| Gaussian curved refraction width | 2.776e-17 | 2e-12 | True |
| Homogeneous absorption and phase | 0 | 2e-12 | True |
| Gaussian beamlet spatial convergence | 2.629e-09 | 2e-06 | True |
| Gaussian beamlet error reduction | 1.273e-09 | 0.0001 | True |
