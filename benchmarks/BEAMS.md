# Beam and dielectric benchmarks

Run `python benchmarks/beam_transport.py`. The script fails if any analytical
check misses its tolerance. It writes CSVs, `summary.json`, `REPORT.md`, and
`convergence.png` to `benchmarks/results/beams`. Values are SI; the phasor
convention is exp(-i omega t). References below are implemented independently
of the solver's Fresnel and curvature transport routines.

| Data file | Independent reference and purpose |
| --- | --- |
| `gaussian_field.csv` | E = E0 exp(ikz)/(1+iz/zR) exp[ik r²/(2(z−izR))], zR = kw0²/2. Samples before, at and after the waist, on and off axis. Checks complex phase as well as amplitude. |
| `gaussian_power_convergence.csv` | Integrate 2πr |E|²/(2η) to radius 5w. Refine 9→513 radial samples. Exact total power is 1 W; omitted tail is exp(−50). |
| `slab_convergence.csv` | Normal-incidence Fabry–Perot sum: t12 t21 p / (1−r21²p²), p = exp(ik0 n d). Includes lossless and complex-index slabs, coherent internal reflections and increasing interaction depth. |
| `branch_pruning_convergence.csv` | Same slab with power cutoff 1e−2→0 at fixed depth. Quantifies field error and discarded branch budget independently of depth. |
| `fresnel_angles.csv` | Textbook real-index s/p coefficients, Snell angle and R+T=1, including Brewster incidence. |
| `concave_mirror.csv` | Spherical mirror q_out^-1 = q_in^-1−2/R, followed by q(z)=q_out+z. Checks width and complex field across focus. |
| `astigmatic_mirror.csv` | Oblique spherical mirror curvature increments −2/(R cos θ) and −2 cos θ/R in tangential/sagittal directions. |
| `absorption.csv` | Homogeneous plane-wave E(z)=exp(ik0 n z), testing attenuation and phase for positive Im(n). |
| `beamlet_convergence.csv` | Convolution decomposition of a Gaussian into narrower Gaussians. Refine 5→25 beamlets per axis; compare off-axis downstream complex fields with the analytical Gaussian, not with fitted launch samples. |

The source decomposition uses the Gaussian convolution identity. For an original
waist w and constituent width b<w, source centers are weighted by
exp[−r0²/(w²−b²)]. Tensor-product trapezoidal quadrature converges to the original
complex beam. Constituent powers cannot be added to measure the power of the
coherent total field. The power quadrature benchmark uses the summed field.

The propagation and local interface maps are analytical paraxial maps; they have
no longitudinal step-size error to converge. Their checks measure agreement
with theory. Actual convergence studies vary branch depth, pruning, receiver
quadrature and beamlet sampling. Low errors establish agreement within this
model, not validity of geometrical optics at arbitrary wavelengths, strong loss,
high numerical aperture, material junctions or clipped beam edges.

Theory references:

- [MIT, Gaussian beams and resonators](https://www.ocw.mit.edu/courses/6-974-fundamentals-of-photonics-quantum-electronics-spring-2006/e9852c138493233bc2813f683da5b199_gaussian_bem_res.pdf)
- [COMSOL, refraction in absorbing media](https://doc.comsol.com/6.3/doc/com.comsol.help.roptics/roptics_ug_optics.6.67.html)

The older C++ spreading benchmark remains a reference for the internal
point-source Jacobian implementation; it is separate from this beam suite.
