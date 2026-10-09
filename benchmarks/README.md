# Analytical spreading benchmarks

`SpreadingBenchmarks.cpp` validates the numerical spreading calculator against
independent closed-form predictions and writes CSV data. `plot_results.py` creates
the comparison figures using Matplotlib in the existing `.venv`.

## Build and run

Configure `cmake-build-tests` using the platform-specific instructions in the
[main README](../README.md#tests): MSVC/vcpkg on Windows, Apple Clang/Homebrew
on macOS, or GCC/distribution packages on Linux. Activate `.venv` first.

Windows Developer PowerShell:

```powershell
./.venv/Scripts/cmake.exe --build cmake-build-tests --target SpreadingBenchmarks
./cmake-build-tests/SpreadingBenchmarks.exe benchmarks/results
./.venv/Scripts/python.exe benchmarks/plot_results.py benchmarks/results
./.venv/Scripts/ctest.exe --test-dir cmake-build-tests --output-on-failure
```

macOS/Linux:

```bash
cmake --build cmake-build-tests --target SpreadingBenchmarks
./cmake-build-tests/SpreadingBenchmarks benchmarks/results
python benchmarks/plot_results.py benchmarks/results
ctest --test-dir cmake-build-tests --output-on-failure
```

CTest writes numerical data into the build directory. Set `BUILD_BENCHMARKS=OFF`
to omit the benchmark executable. Plotting is a separate Python operation.

## Definition and normalization

The calculator perturbs a unit launch direction in an orthonormal tangent
basis, replays its fixed reflection sequence, and intersects neighbouring rays
with a fixed plane through the receiver normal to the central arrival direction.
The two receiver transverse coordinates define the 2-by-2 angular Jacobian J.

- Area per source solid angle: `A = |det J|` (square metres per steradian locally).
- Electric field multiplier: `g = r_ref / sqrt(A)`.
- Default source reference distance: `r_ref = 1 m`.

The source field must be interpreted as the homogeneous outgoing point-source
field at `r_ref`, along the launch direction. Reflection polarization and optical
phase are separate calculations. Applying this factor once accounts for the
entire path's geometric spreading; multiplying an extra `1/L` would double-count it.

The implementation compares central differences at h and h/2 and halves h when
the neighbourhood is invalid or derivatives disagree. It checks the matrix and
both singular values. Results include the actual step, derivative difference,
tangent bases and diagnostic status. A central ray must reach the point receiver
within tolerance. Near-zero singular values produce `Caustic` without a field factor.

## Independent analytical predictions

Let a be distance from source to reflection, b distance from reflection to
receiver, R the positive geometric radius, and theta incidence angle from the
normal. The predicted singular values are the absolute values of the following
spreading lengths, sorted in descending order. Their product is A.

| Case                                                       | First spreading length | Second spreading length |
|------------------------------------------------------------|------------------------|-------------------------|
| Direct ray at distance L                                   | L                      | L                       |
| Any valid planar reflection sequence, total length L       | L                      | L                       |
| Convex sphere, normal incidence                            | a+b+2ab/R              | a+b+2ab/R               |
| Convex cylinder, normal incidence                          | a+b+2ab/R              | a+b                     |
| Convex sphere, oblique incidence                           | a+b+2ab/(R cos theta)  | a+b+2ab cos theta/R     |
| Convex cylinder, incidence plane perpendicular to its axis | a+b+2ab/(R cos theta)  | a+b                     |
| Concave sphere, normal incidence                           | a+b-2ab/R              | a+b-2ab/R               |

A two-bounce convex-sphere case additionally validates accumulated curvature.
Its independent analytical transfer matrix is
`T = P(b) M(R2) P(d) M(R1) P(a)`, where
`P(l) = [[1,l],[0,1]]` and `M(R) = [[1,0],[2/R,1]]`.
Both spreading lengths are `T[0,1]`. The scene uses a=1 m, d=3 m,
R1=1 m and R2=0.7 m; the receiver stays between the reflectors.

These are differential ray results: the finite angular bundle has higher-order
aberrations, but its derivative at the central ray has the stated analytical value.
For normal incidence, the reflected angular slope is `1 +/- 2a/R`, giving
`B = a + b(1 +/- 2a/R)`. The plus sign is for exterior convex reflection; the minus
sign is for interior concave reflection. Oblique spherical reflection has distinct
tangential and sagittal wavefront curvatures.

References:

- [Bremmer, Reflection from an arbitrarily curved surface (1982)](https://doi.org/10.1029/RS017i005p01117).
- [NASA reflector analysis, equations 2.31–2.32](https://ntrs.nasa.gov/api/citations/19880016378/downloads/19880016378.pdf),
  expressing field divergence through the reflected wavefront's two principal radii.

The sweeps use a=2 m and R=1 m for convex reflections, with 45-degree and 70-degree
spherical cases and a 45-degree cylindrical case. The concave sphere uses a=3 m,
R=2 m, with its focus at b=1.5 m. At the focus the plotted curves have a gap and the
calculator is separately required to report a caustic. A concave cylinder also
tests a one-axis caustic. Every regular sweep point is repeated after a rigid
rotation and translation.

Additional checks cover source normalization, a receiver miss, a central blocker,
an aperture edge, invalid options, unstable derivatives, and reducing the angular
step around a blocked neighbouring ray. A separate oblique-sphere step sweep
requires approximately second-order convergence of the area error.

## Outputs and scope

`results/` contains ten comparison PNGs, a convergence PNG, `spreading.csv`,
`convergence.csv`, and `summary.txt`. Each comparison shows analytical field factor
and Jacobian samples alongside the maximum relative error across both singular
values, area and field factor. Numerical checks use a relative tolerance of 2e-6.
Error panels label their scientific scaling explicitly; CSV errors are unscaled.
Plotting floors zero errors at 1e-16 to keep the logarithmic axis defined. The
summary's maximum also includes the rotated and translated scenes.

The concave sweep validates amplitude magnitude on both sides of a focus. The
calculator does not track earlier caustic crossings or their phase shifts, and
does not regularize the divergent geometrical-optics amplitude at a focus.
