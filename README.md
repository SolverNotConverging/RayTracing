# Ray tracing solver

The C++ API launches rays, refines received reflection sequences to a point
receiver, deduplicates paths, calculates Jacobian spreading, and reconstructs
complex electric fields and a coherent impulse response. Geometry uses specular
PEC reflections in a homogeneous, nondispersive medium.

## Library and application

The root CMake project builds `app/main.cpp`. All solver implementation and public
headers live in `RayTracing/`, with its own CMake fragment. Three static libraries
expose these targets:

| Target                      | API                                                        |
|-----------------------------|------------------------------------------------------------|
| `RayTracing::RayTracing`    | Antennas, scene, tracing, refinement, spreading and fields |
| `RayTracing::IO`            | Core plus HDF5 and impulse CSV persistence                 |
| `RayTracing::Visualization` | IO plus VTK views and Matplot++ plots                      |

The executable target is `RayTracingApp`, with output name `RayTracing`.
Configuration is ordinary C++ in `main.cpp`; `trace_options.json` was removed.

```cpp
#include <RayTracing/Solver.hpp>
#include <RayTracing/ResultIO.hpp>
#include <RayTracing/Visualization.hpp>

rt::SolverConfig config;
config.frequencyHz = 77e9;
config.rayCount = 2000;
config.maxReflections = 8;
config.maxDistance = 20;

rt::Scene scene;
scene.add(rt::Sphere{{0.7, -1.2, 0}, 0.3});
rt::Transmitter tx{{0, 0, 0}, rt::ShortDipole{1e-4, 1.0}};
rt::Receiver rx{{-0.5, 1.4, 0}, rt::Isotropic{}};

auto result = rt::solve(scene, tx, rx, config);
rt::save_h5(result, "results/simulation.h5");
rt::save_impulse_csv(result.impulseResponse, "results/impulse_response.csv");
rt::plot(result.impulseResponse, "results/impulse_response.png", false);
rt::visualize(result);  // rt::visualize(result, 0) selects one path.
```

`result.rays` contains converged distinct paths, spreading and optional fields.
`result.candidates` retains refinement diagnostics and coarse received paths.
Invalid spreading leaves the path available for inspection with no field.
Geometry discovery continues through source-pattern nulls. An additional exact
direct launch is included, so `launchedRays = rayCount + 1`; it undergoes the same
visibility checks and deduplication as sampled launches.

## Tx and Rx antennas

Both endpoints use `rt::Antenna`. Its `orientation` maps local antenna coordinates
into world coordinates and must be a proper rotation. Positions are separate.

| Model                 | Parameters and local convention                                                             |
|-----------------------|---------------------------------------------------------------------------------------------|
| `Isotropic`           | Spherical theta/phi or circular polarization; complex source amplitude                      |
| `ShortDipole`         | Effective uniform-current length and complex current; local z                               |
| `ThinWireDipole`      | Length and complex feed current; sinusoidal current, centre fed, local z                    |
| `RectangularAperture` | Width, height and uniform complex tangential Ex/Ey; local +z radiation, infinite PEC baffle |
| Imported pattern      | CST `.ffs` or HFSS `.ffd`, phase convention and field scaling; power reference for Rx       |

An isotropic Rx still selects polarization, with unit directional amplitude.
Its source amplitude setting applies to Tx use. Its spherical polarization basis
has a pole convention. Circular labels mean `(e_theta +/- i e_phi)/sqrt(2)` under
the solver's time convention. `receiveCalibration` supplies an explicit complex
Rx port amplitude/phase calibration.

```cpp
rt::FarfieldImportOptions options;
options.inputConvention = rt::PhasorConvention::PositiveTime;
tx.antenna = rt::load_farfield("simple_patch.ffs", options);
rx.antenna = rt::load_farfield("simple_patch.ffs", options);
rx.antenna.orientation =
    Eigen::AngleAxisd(0.5, rt::Vec3::UnitZ()).toRotationMatrix();
```

The supplied CST version 3.0 file has five frequencies from 76 to 78 GHz, with
361 azimuth and 181 elevation samples each. At 77 GHz it records 0.3823696 W
radiated, 0.4978351 W accepted and 0.5 W stimulated power. Both complex spherical
field components, powers and exported coordinate metadata are preserved. Exported
axes initialize antenna orientation; exported position is metadata, while the
endpoint position controls placement.

CST version 3.0 Farfield Source and rectangular theta/phi HFSS exports are
supported. HFSS accepts frequency-independent data or `Frequencies`/`Frequency`
blocks in the [documented FFD layouts](https://help.agi.com/stk/12.8.0/Content/comm/complexANSYSffdSamples.htm).
HFSS lacks input-power metadata: supply `inputPowerWatts` for Rx use. Select
`receivePowerReference` as accepted, radiated or stimulated power (default:
accepted). Supplied HFSS input power populates accepted/stimulated references;
a radiated reference requires separately supplied pattern metadata.

Fields are canonical Cartesian complex `rE` coefficients in volts, with outgoing
radial propagation removed. `coefficientScale` converts export units/reference
radius to that convention. Positive-time input is conjugated into
`exp(-i omega t)`; choose the option matching the export, since files do not
declare it. Angular interpolation is Cartesian bilinear followed by transverse
projection; frequency interpolation is complex linear. Complete azimuth grids
wrap, and coverage is enforced without extrapolation. Custom sampled analytical
patterns can be supplied through `FarfieldData` with explicit grids, Cartesian
coefficients and powers.

VTK replaces each non-isotropic endpoint sphere with its oriented pattern.
Its radial display uses `|rE| / max |rE|`; `ViewOptions::patternScale` sets the peak
display radius in metres and leaves the physical fields unchanged.

## Fields and receive response

For path p with physical length L, medium index n, angular Jacobian J and PEC
normals n_j, the solver uses

$$P_j=-I+2\mathbf n_j\mathbf n_j^T,\qquad
\mathbf E_p=\frac{P_N\cdots P_1\mathbf F_t(\hat{\mathbf d}_t,f)}
{\sqrt{|\det J_p|}}\exp(+i\,2\pi f nL_p/c).$$

F_t is the source `rE` vector. The low-level reference-distance factor cancels
its source-field conversion. One common `config.commonSourceReference` normalizes
the fields, preserving directional amplitudes. The medium uses
`n=sqrt(epsilon_r mu_r)` and `eta=376.730313668 sqrt(mu_r/epsilon_r)` ohms.
The existing `c=3e8 m/s` is retained.

For a patterned reciprocal Rx, the outward look direction is the negative ray
arrival direction, and

$$\mathbf w_r=\sqrt{\frac{4\pi}{2\eta P_{\rm ref}}}\mathbf F_r,\qquad
a_p=\gamma_r\mathbf w_r^T\mathbf E_p/E_{\rm common}.$$

This bilinear reciprocity contraction uses the transmit-pattern representation
and preserves complex phase and polarization. A unit isotropic polarization
vector replaces w_r for an isotropic Rx. The low-level fixed-analyzer impulse API
retains its Hermitian analyzer convention.

The output is a normalized field channel with reciprocal Rx gain weighting.
Absolute terminal voltage, mismatch and S-parameters require a specified port or
effective-length calibration. Incident vector fields are retained separately.
`rt::update_receiver(result, antenna)` reapplies an Rx pattern/orientation without
tracing again.

$$h(t)=\sum_p a_p\delta(t-\tau_p),\quad \tau_p=nL_p/c,\qquad
H(\nu)=\sum_p a_p\exp(+i\,2\pi\nu\tau_p).$$

Carrier phase is already in a_p. Delays within `delayToleranceSeconds` of each
group's earliest delay merge coherently. The default 1e-13 s tolerance handles
numerical equality; bandwidth and pulse shaping are separate inputs. Plots show
absolute delay, linear magnitude and wrapped phase. Phase below 1e-12 times the
strongest tap is omitted. Receiver caustics are detected; earlier caustic crossings
and their phase shifts remain outside the current model.

## Save, reload and view

```cpp
auto saved = rt::load_h5("results/simulation.h5");
rt::visualize(saved);
rt::visualize_h5("results/simulation.h5");
auto response = rt::load_impulse_csv("results/impulse_response.csv");
rt::plot(response);
rt::plot_csv("results/impulse_response.csv", "replotted.png", false);
```

HDF5 schema version 1 stores settings, all geometry, Tx/Rx definitions and embedded
patterns, candidate/refined paths, spreading diagnostics, vector fields, coherent
taps and tracing counts. Groups are `/meta`, `/settings`, `/scene`, `/transmitter`,
`/receiver`, `/paths`, `/candidates`, `/response` and `/diagnostics`. Complete
records and definitions are UTF-8 JSON datasets. Native numeric datasets expose
pattern coefficients `(frequency, phi, theta, xyz, real/imag)`, path vertices and
offsets, lengths, complex fields, Jacobians and taps. Reloading needs no original
antenna file.

CSV preserves carrier frequency, phase convention, Rx model, complex taps, vector
sums and contributing path indices. Magnitude/phase columns aid inspection;
undefined phases are blank. Path indices map through `result.responseRayIndices`
into `result.rays`.

## Build and run

Dependencies are Eigen3, nlohmann_json, HDF5, VTK and Matplot++ (vcpkg packages
`eigen3`, `nlohmann-json`, `hdf5`, `vtk`, `matplotplusplus`). Gnuplot must be on
`PATH` for plotting. Use a compiler developer shell and the existing CMake/vcpkg
configuration.

```powershell
cmake --build cmake-build-debug --target RayTracingApp
./cmake-build-debug/RayTracing.exe
./cmake-build-debug/RayTracing.exe --no-gui --output results --screenshot results/scene.png
./cmake-build-debug/RayTracing.exe --isotropic --no-gui --output results/isotropic
./cmake-build-debug/RayTracing.exe --load results/simulation.h5 --no-gui --output results/reloaded
./cmake-build-debug/RayTracing.exe --plot-csv results/impulse_response.csv --no-gui --output results/replotted
```

The demo uses the supplied CST Tx and an isotropic Rx. It exports `simulation.h5`,
`impulse_response.csv` and a magnitude/phase PNG after tracing. Interactive mode
opens Matplot++ first; close it to proceed to VTK. `--no-gui` still exports the
plot, and `--screenshot` exports VTK without an interactive window in that mode.
For a numerical build without VTK/Matplot++, configure
`-DRAYTRACING_BUILD_VISUALIZATION=OFF -DBUILD_BENCHMARKS=OFF`.

## Validation

```powershell
cmake --build cmake-build-debug
ctest --test-dir cmake-build-debug --output-on-failure
```

Six suites cover geometry, sequence replay, refinement, field/impulse math,
antennas/imports/API/persistence and analytical spreading. Antenna tests use the
real CST file and synthetic independent/multifrequency HFSS fixtures. They check
amplitude/phase, polarization, rotations, nulls, shared source normalization,
Rx replacement and nonempty/empty/zero-field HDF5/CSV round trips. Analytical
spreading formulas and comparison plots are in [benchmarks/README.md](benchmarks/README.md).
