# RayTracing

Python simulations backed by a C++20 numerical solver. The solver launches rays,
refines specular PEC paths to a point receiver, deduplicates them, calculates
Jacobian spreading, and reconstructs complex electric fields and coherent impulse
responses in a homogeneous, nondispersive medium. Matplotlib plots responses;
PyVista displays geometry, antenna patterns, paths and polarization.

## Install in the existing environment

The native build requires MSVC, Eigen3, nlohmann_json and HDF5. This repository's
Windows setup uses `C:/opt/vcpkg`. Use a Visual Studio developer PowerShell and
the existing `.venv` (Python 3.14). All Python dependencies stay in that environment:

```powershell
./.venv/Scripts/python.exe -m pip install "scikit-build-core>=0.11" "pybind11>=3" "cmake>=4.2" ninja build
$env:CMAKE_ARGS = '-DCMAKE_TOOLCHAIN_FILE=C:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake'
./.venv/Scripts/python.exe -m pip install --no-build-isolation ".[dev,notebooks]"
./.venv/Scripts/python.exe -m ipykernel install --sys-prefix --name raytracing --display-name "RayTracing (.venv)"
./.venv/Scripts/python.exe -m jupyterlab
```

Open [example1.ipynb](example1.ipynb) or [example2.ipynb](example2.ipynb) with the
**RayTracing (.venv)** kernel. Example 1 contains the original mixed geometry and
horizontal CST transmitter. Example 2 has a cylinder of radius 0.127 m from z=0
to -0.5 m, open at the top, with a bottom disk. Both CST patches face downward
at (0.005, 0, 0) m. The common 5 mm shift avoids the exact axial caustic while
retaining monostatic reception. Set `source_shift = 0` to inspect that singular
case. Each notebook saves its HDF5, CSV and images under `results/exampleN`.

There is one complete Python package. Matplotlib and PyVista are installed as
Python dependencies and imported only when rendering. PyVista supplies its own
Python VTK dependency; the C++ solver neither includes nor links VTK. The Windows
wheel bundles its HDF5 runtime and CST example pattern. The C++ application and
native rendering interfaces have been removed in this breaking release.

## Python API

```python
import raytracing as rt

scene = rt.Scene()
scene.add(rt.Sphere([0.7, -1.2, 0], 0.3))
patch = rt.load_farfield(rt.example_pattern())
patch.rotate([0, 1, 0], -90)  # Native -z beam rotated toward +x.
tx = rt.Transmitter([0, 0, 0], patch)
rx = rt.Receiver([-0.5, 1.4, 0], rt.Isotropic())
config = rt.SolverConfig(frequency_hz=77e9, ray_count=2000, max_reflections=8)
result = rt.solve(scene, tx, rx, config)

vertices = result.rays[0].vertices  # N-by-3 NumPy coordinates.
response = rt.frequency_response(result.impulse_response)
rt.save_h5(result, "simulation.h5")
rt.save_impulse_csv(result.impulse_response, "impulse_response.csv")
figure = rt.plot(result.impulse_response, "impulse_response.png")
viewer = rt.visualize(result)
viewer.screenshot("scene.png")
viewer.close()
```

Coordinates and lengths use metres; frequencies use Hz; rotations use degrees.
Vector and matrix results are NumPy arrays, fields are complex-valued. Configurations
and antenna definitions are editable; solve results expose numerical diagnostics
through read-only properties. `solve` releases the Python GIL while computing.

| API | Purpose |
| --- | --- |
| `Scene.add(Sphere / Rectangle / Disk / Cylinder / Triangle)` | Construct geometry |
| `Isotropic`, `ShortDipole`, `ThinWireDipole`, `RectangularAperture` | Analytical antenna models |
| `load_farfield(path, options)` | Import CST `.ffs` or rectangular-grid HFSS `.ffd` |
| `solve(scene, tx, rx, config)` | Compute paths, fields and impulse response |
| `evaluate_sequence`, `refine_path`, `calculate_spreading`, `deduplicate_paths` | Lower-level numerical operations |
| `update_receiver(result, antenna)` | Apply a new Rx antenna without retracing |
| `save_h5`, `load_h5`, `save_impulse_csv`, `load_impulse_csv` | Persist and reload results |
| `frequency_response(response, offset_hz=0)` | Evaluate the complex channel response |
| `plot`, `plot_csv` | Return a Matplotlib Figure; optionally save it |
| `visualize`, `visualize_h5` | Return a PyVista Plotter |

`result.rays` contains distinct converged paths with geometry, spreading and optional
fields. `result.candidates` retains coarse paths and refinement failures. Invalid
spreading retains the geometry without assigning a field. The reception sphere
collects candidates along every ray segment; it does not absorb the ray at the
first near-pass. Distinct endpoints get one additional exact direct launch.
Colocated endpoints receive reflected returns without a zero-length direct ray.
Refinement defaults to a 1e-10 m receiver tolerance. Cylinder rim hits are marked
`UNRESOLVED_EDGE`; receiver caustics have `CAUSTIC` spreading status.

## Antennas and visualization

The supplied CST patch contains five frequencies from 76 to 78 GHz. Its main
beam points along -z. Antenna orientation maps local coordinates to world
coordinates; `rotate(axis, degrees)` composes a right-hand world-axis rotation.
Endpoints copy antenna definitions, so their orientations can be changed independently.
Imported pattern samples are shared and read-only.

`FarfieldImportOptions` selects `input_convention`, `receive_power_reference`,
`input_power_watts`, `coefficient_scale` and `honor_export_axes`. The default CST
input convention is positive-time, conjugated into the solver's `exp(-i omega t)`
convention. HFSS exports need an explicit input power for Rx use. Fields are
Cartesian complex `rE` coefficients in volts. Angular interpolation is Cartesian
bilinear followed by transverse projection; frequency interpolation is complex
linear. Complete azimuth grids wrap; partial coverage does not extrapolate.

`rt.plot` returns an ordinary Matplotlib Figure; use its axes to customize it and
`matplotlib.pyplot.close(figure)` when finished. PyVista views accept
`options=rt.ViewOptions(...)`, including `ray_indices`, `pattern_scale`,
`field_dynamic_range_db`, `show_polarization` and `show_unresolved_candidates`.
The default off-screen Plotter supports notebook display and image export:

```python
saved = rt.load_h5("simulation.h5")
viewer = rt.visualize(saved, options=rt.ViewOptions(ray_indices=[0]))
viewer.show(jupyter_backend="static")
```

You can also use PyVista's interactive desktop view with
`ViewOptions(off_screen=False, notebook=False)` and `viewer.show()`. This is a
runtime presentation choice within the same package.

Patterns use normalized radial magnitude for display, preserving world orientation.
Red denotes Tx and green Rx. Paths are colored by incident Rx field magnitude,
`20 log10(|E_p| / max |E|)`, before receive weighting or coherent addition. The
reference is the strongest full-result path even when displaying a subset.
Unavailable fields appear grey. Polarization ellipses have bounded display size.
Unresolved corner and edge candidates can be displayed in orange.

## Fields and receive response


For path p with physical length L, medium index n, angular Jacobian J and PEC
normals n_j, the solver uses

$$P_j=-I+2\mathbf n_j\mathbf n_j^T,\qquad
\mathbf E_p=\frac{P_N\cdots P_1\mathbf F_t(\hat{\mathbf d}_t,f)}
{\sqrt{|\det J_p|}}\exp(+i\,2\pi f nL_p/c).$$

F_t is the source `rE` vector. The low-level reference-distance factor cancels
its source-field conversion. One common `config.common_source_reference` normalizes
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
`rt.update_receiver(result, antenna)` reapplies an Rx pattern/orientation without
tracing again.

$$h(t)=\sum_p a_p\delta(t-\tau_p),\quad \tau_p=nL_p/c,\qquad
H(\nu)=\sum_p a_p\exp(+i\,2\pi\nu\tau_p).$$

Carrier phase is already in a_p. Delays within `delay_tolerance_seconds` of each
group's earliest delay merge coherently. The default 1e-13 s tolerance handles
numerical equality; bandwidth and pulse shaping are separate inputs. Plots show
absolute delay, linear magnitude and wrapped phase. Phase below 1e-12 times the
strongest tap is omitted. Receiver caustics are detected; earlier caustic crossings
and their phase shifts remain outside the current model.

## Persistence

HDF5 schema version 1 stores settings, geometry, both antenna definitions and
embedded patterns, candidate/refined paths, spreading diagnostics, vector fields,
coherent taps and tracing counts. Reloading requires no original antenna file.
CSV preserves carrier frequency, phase convention, Rx model, complex taps, vector
sums and contributing path indices. These indices map through
`result.response_ray_indices` into `result.rays`.

## Build and validation

CMake exposes `RayTracing::RayTracing` for numerics and `RayTracing::IO` for
persistence. Python bindings are in `python/bindings`; Python rendering lives in
`python/raytracing/visualization.py`. There are no native visualization targets.

```powershell
./.venv/Scripts/python.exe -m build --wheel --no-isolation
./.venv/Scripts/cmake.exe -S . -B cmake-build-tests -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=C:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake
./.venv/Scripts/cmake.exe --build cmake-build-tests
./.venv/Scripts/ctest.exe --test-dir cmake-build-tests --output-on-failure
./.venv/Scripts/python.exe -m pytest tests/python
```

Seven C++ suites cover surfaces, sequence replay, refinement, fields, antennas and
persistence, cylinder returns and analytical spreading. Python tests exercise the
bindings, NumPy ownership, pattern rotation, shifted cylinder, persistence, plots
and installed-package runtime. Both notebooks are executable integration examples.
See [benchmarks/README.md](benchmarks/README.md) for analytical formulas and
Matplotlib comparison plots.
