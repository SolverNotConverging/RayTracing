# RayTracing

A standalone C++20 point-ray and Gaussian-beam numerical library with a Python API, Matplotlib
impulse plots and PyVista ray visualization. The solver discovers and refines
specular paths, calculates spreading and polarized fields, and returns coherent
receiver impulse taps. Graphics libraries belong to applications, outside the
numerical library.

Both solvers run in the same C++ core. Python provides data adapters, persistence
and visualization; it calls the native numerical implementation for both APIs.

## Standalone C++ build

Requires CMake 4.2+, a C++20 compiler and Eigen3. Python, HDF5, JSON, Matplot++
and VTK are not required for a core build. In a Windows Developer PowerShell,
install Eigen with vcpkg and pass your toolchain path:

```powershell
cmake -S . -B build/native -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=C:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build/native
ctest --test-dir build/native --output-on-failure
./build/native/examples/cpp/point_ray_example.exe --rays 20000
cmake --install build/native --prefix ./build/sdk --component Development
```

On Linux/macOS use the same CMake commands with your Eigen installation, omit the
vcpkg toolchain if unnecessary, and omit `.exe`. The core is a static library;
Eigen and C++20 are public requirements. Build consumers with compatible compiler,
Eigen settings and C++ runtime. No cross-compiler binary ABI is promised.

To consume an installation:

```cmake
cmake_minimum_required(VERSION 4.2)
project(MySimulation LANGUAGES CXX)
find_package(RayTracing 0.5 CONFIG REQUIRED COMPONENTS Core)
add_executable(simulation main.cpp)
target_link_libraries(simulation PRIVATE RayTracing::Core)
```

Set `CMAKE_PREFIX_PATH` to the SDK installation prefix. Alternatively,
`add_subdirectory(path/to/RayTracing)` provides the same target without installing.
Examples and tests are off by default when embedded in another project.

```cpp
#include <RayTracing/Solver.hpp>
#include <cmath>

int main() {
    rt::Scene scene;
    rt::Transmitter tx{Vec3::Zero(), rt::Isotropic{}};
    rt::Receiver rx{{3, 0, 0}, rt::Isotropic{}};
    const auto result = rt::solve(scene, tx, rx);
    const auto channel = rt::frequency_response(result.impulseResponse);
    return result.rays.empty() || !std::isfinite(std::abs(channel));
}
```

Imported antenna patterns are loaded explicitly at runtime using `rt::load_farfield(path)`
in C++ or `rt.load_farfield(path)` in Python. No pattern data is compiled into the
library or bundled into the wheel. CST/HFSS files in `antenna_patterns/` are runtime
test fixtures. Examples use analytical antennas.

Use SI units. The field convention is `exp(-i omega t)`. Results own their scene,
settings and antenna snapshots. `responseRayIndices` maps impulse path indices to
solved rays; a refined geometric path can have no field when spreading is singular.
Invalid inputs throw exceptions, while numerical path statuses describe unsuccessful
or singular candidates. Low-level geometry routines require validated geometry.

## C++ persistence and graphics

| CMake option | Default | Effect |
| --- | --- | --- |
| `RAYTRACING_BUILD_PYTHON` | OFF | Full Python API; requires IO |
| `RAYTRACING_BUILD_IO` | Same as Python option on first configuration | HDF5/CSV library |
| `RAYTRACING_BUILD_EXAMPLES` | ON for standalone builds | C++ room example |
| `RAYTRACING_EXAMPLE_MATPLOT` | OFF | Matplot++ impulse plotting in example |
| `RAYTRACING_EXAMPLE_VTK` | OFF | Native VTK ray viewer in example |
| `RAYTRACING_BUILD_BENCHMARKS` | OFF | Analytical spreading benchmarks |
| `BUILD_TESTING` | ON for standalone builds | Numerical tests; IO checks when enabled |

Enable optional native features explicitly:

```powershell
cmake -S . -B build/graphics -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=C:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake -DRAYTRACING_BUILD_IO=ON -DRAYTRACING_EXAMPLE_MATPLOT=ON -DRAYTRACING_EXAMPLE_VTK=ON
cmake --build build/graphics
./build/graphics/examples/cpp/point_ray_example.exe --output build/example-output --plot build/example-output/impulse.png --screenshot build/example-output/scene.png
./build/graphics/examples/cpp/point_ray_example.exe --view
```

Install `nlohmann-json` and `hdf5` for IO, `matplotplusplus` for plotting, and `vtk`
for the native viewer. The Matplot++ default backend also needs gnuplot at runtime.
These native graphics dependencies may require building through your package
manager. Python visualization uses the installed Matplotlib/PyVista packages and
does not require enabling either C++ graphics option.

The C++ example uses the same mixed-geometry room as `example1.py`, with an
analytical half-wave dipole transmitter and isotropic receiver. `--plot` saves magnitude and phase;
`--screenshot` renders offscreen unless `--view` is also supplied. Cyan paths have
fields, orange paths have invalid spreading, Tx is red and Rx is green. Use
`--help` for all options. No files or windows are created without those flags.

The example directory is also a standalone SDK consumer:

```powershell
cmake -S examples/cpp -B build/consumer -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$PWD/build/sdk" -DCMAKE_TOOLCHAIN_FILE=C:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build/consumer
```

For persistence, request `find_package(RayTracing CONFIG REQUIRED COMPONENTS IO)`,
link `RayTracing::IO`, and include `<RayTracing/ResultIO.hpp>`. The consumer must
enable C as well as C++ because CMake's HDF5 finder checks its C library. An
installed package loads IO dependencies only when the IO component is requested.

## Python installation

Requires CPython 3.14+, the C++ build tools above, Eigen3, nlohmann-json and HDF5.
The Python package always includes NumPy, Matplotlib and PyVista dependencies.
Native VTK and Matplot++ are not used by the Python build.

```powershell
# In an x64 Visual Studio Developer PowerShell, with your environment active:
$env:CMAKE_GENERATOR = 'Ninja'
$env:CMAKE_ARGS = '-DCMAKE_TOOLCHAIN_FILE=C:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows -DCMAKE_CXX_SCAN_FOR_MODULES=OFF'
python -m pip install .
```

The wheel build explicitly enables Python and IO and installs only the Python
component, excluding C++ SDK headers and archives. For development, install
`.[dev]`. The local `.\install.ps1` shortcut remains usable with the local preset.
On Linux/macOS install the corresponding native libraries and configure
`CMAKE_PREFIX_PATH` as needed.

The point-ray API is `from raytracing import point_ray as rt`. The top-level
`raytracing` API continues to expose Gaussian-beam transport. PyCharm can mark
`python/` as a Sources Root; source imports find the native extension in the
matching installed wheel. Rebuild/reinstall after C++ changes.

## Main examples: multipath impulse responses

| Example | Scene | Launch rays |
| --- | --- | --- |
| [example1.py](example1.py) / [notebook](example1.ipynb) | Mixed-geometry room; analytical half-wave dipole along z; isotropic receiver | 20,000 |
| [example2.py](example2.py) / [notebook](example2.ipynb) | Open PEC cylinder and bottom disk; colocated horizontal analytical dipoles offset 6 mm from the axis | 200,000 |

Both scenes use explicit PEC surface materials, a 77 GHz carrier, up to eight
reflections and a 20 m path limit. Run from the repository with its `.venv`:

```powershell
.\.venv\Scripts\python.exe example1.py
.\.venv\Scripts\python.exe example2.py
```

Each script opens an interactive scene. Add `--no-show` to generate results off
screen. The notebooks display the impulse plots and scene screenshots inline;
select the repository's `.venv` kernel.

Use `--rays N` to change the discovery budget and `--output DIR` to choose the
output directory. By default, each example writes these files under `results/example1` or `results/example2`:

- `impulse_response.png`: coherent tap magnitude and wrapped phase versus absolute delay.
- `impulse_response.csv`: complex impulse taps and their delays.
- `simulation.h5`: saved simulation snapshot.
- `scene.png`: antenna, geometry and path visualization.

The point-ray HDF5 format stores large path and discovery-candidate records in
chunked byte datasets. These records can exceed 64 MiB and round-trip through
`save_h5` and `load_h5`; the reader also accepts legacy scalar path records.
Other JSON metadata fields retain a 64 MiB limit, enforced by both writer and reader.

The point-ray API uses `Scene.add(shape, material)` and
`solve(scene, transmitter, receiver, config)`. To replot a saved response:

```python
from raytracing import point_ray as rt

result = rt.load_h5("results/example1/simulation.h5")
figure = rt.plot(result.impulse_response)
figure.show()
channel_at_carrier = rt.frequency_response(result.impulse_response)
```

The taps retain coherent carrier phase. `frequency_response(response, offset_hz)`
evaluates the envelope response at a frequency offset from that carrier using the
computed taps. Launch-ray density controls path discovery; refine it when checking
convergence. Singular receiver caustics have no valid point-ray field contribution.

## Native Gaussian-beam API

`<RayTracing/BeamSolver.hpp>` exposes the Gaussian solver in `rt::beam`, linked
through the same `RayTracing::Core` target. It includes material evaluation,
volume insertion priority, termination, Gaussian/plane-wave/aperture launch,
beamlet decomposition, Fresnel splitting, differential curved-interface maps,
propagation, field sampling and branch power diagnostics. It depends only on Eigen
and the C++ standard library. It never imports Python or creates plotting windows.

```cpp
#include <RayTracing/BeamSolver.hpp>

namespace b = rt::beam;
b::Scene scene;
scene.add(b::Box{{0, 0, .5065}, {2, 2, .0065}},
          b::Material{"dielectric", {4, .02}, 1, 0, "glass"}, true);
b::GaussianBeam beam;
beam.waistRadius.setConstant(.1);
b::SolverConfig config;
config.frequencyHz = 10e9;
config.maxDistance = 4;
const auto result = b::solve(scene, b::Source{beam}, config);
const Vec3C field = result.field_at({.03, .01, 1});
```

Run `build/native/examples/cpp/gaussian_beam_example` (with `.exe` on Windows) for
complex-field CSV samples through a lossy dielectric slab. The example also builds
from the standalone installed-package consumer directory. Plot its output in an
application; the numerical example has no rendering dependency.

Use a `std::vector<b::Source>` for coherent source superposition. An
`ApertureSource::field` callback supplies Cartesian complex field samples; C++
constructs the grid, fits coefficients and launches beamlets. The
`ApertureSource::from_antenna(...)` factory samples analytical or runtime-imported
antenna patterns in the vacuum far zone. `solve_states` accepts custom launch
states. Launch states, geometry and settings are copied into numerical results.

Python keeps the existing `Scene`, `GaussianBeam`, `PlaneWave`, `ApertureSource`,
`SolverConfig`, `Simulation`, plotting and JSON APIs. Python callables provide
aperture samples while holding the GIL; the C++ branch loop and field queries
release it. Python result data remain serializable and editable; field queries
convert the current data into a native snapshot, so edits cannot leave a stale
cached numerical result. Batch `fields_at` amortizes this conversion.

All model limits below still apply: this port preserves the paraxial, weak-loss
model and does not add clipping, diffraction or receive-mode overlap. The Gaussian
solver returns field samples and branch trees; impulse taps remain the point-ray
solver's output.

## Gaussian-beam transport

The following slab snippet uses the top-level Gaussian-beam API. Reproducible
slab, mirror and refraction validations are maintained in
`benchmarks/beam_transport.py` rather than duplicate example scripts/notebooks:

```python
import raytracing as rt

scene = rt.Scene()
scene.add_volume(
    rt.Box(center=[0, 0, 0.51], half_lengths=[0.5, 0.5, 0.01]),
    rt.Material.dielectric(4 + 0.02j, name="lossy glass"),
)
scene.add_termination(
    rt.Box([0, 0, 0.8], [2, 2, 1.8]), volume=True, on="exit",
)
beam = rt.GaussianBeam(
    waist_position=[0, 0, 0.8], direction=[0, 0, 1],
    waist_radius=0.04, power=1.0, polarization=[1, 0, 0],
    launch_distance=-0.8,
)
config = rt.SolverConfig(frequency_hz=77e9)
result = rt.solve(scene, beam, config)
print(result.field_at([0.01, 0, 1.0]))  # Cartesian complex E, V/m
print(result.power_balance)
print(result.warnings)
rt.save_result(result, "beam.json")
rt.visualize(result)
```

All dimensions are metres. `waist_radius` is the 1/e field radius (1/e^2 intensity
radius). It accepts one radius or two transverse radii. `launch_distance` is signed
axial distance from the waist; a negative value starts before the focus. An
optional `transverse_axis` defines the first elliptical waist axis. The source
waist specifies the incident beam in its launch medium; interfaces subsequently
move and reshape its focus. Source and launch plane must lie in a homogeneous
lossless launch region. `PlaneWave(position, direction, polarization, amplitude)`
is available for interface/slab calculations.

## Gaussian-beam geometry and materials

| Method | Meaning |
| --- | --- |
| `Scene(background=Material.dielectric(...))` | Homogeneous background |
| `add_volume(shape, material)` | Closed material region |
| `add_surface(shape, material)` | Opaque PEC/metal boundary |
| `add_termination(shape, volume=False)` | Stop every ray/beam axis hitting this surface |
| `add_termination(shape, volume=True, on="entry")` | Stop on entering a volume, including launches inside it |
| `add_termination(shape, volume=True, on="exit")` | Stop on leaving an enclosing region |

Volumes support oriented `Box`, `Sphere` and capped `Cylinder`. Rectangles,
disks, triangles and open cylinders are surfaces. Dielectrics must be volumes.
Both PEC and finite-conductivity metal are opaque: neither produces transmitted
branches, even when supplied as volumes. Metal uses lossy Fresnel reflection;
PEC uses its exact reflection limit. There is no skin-depth propagation or
thin-metal-film transmission model.

The newest material volume wins at an overlapping point. Within the interior of
coincident opaque surfaces, the newest surface wins. Materials are resolved on
both sides of every hit; equal-material internal boundaries are skipped.
A shared smooth face between two volumes is a valid interface. Edges, rims,
corners and competing normals produce an `ambiguous` terminal. Termination
boundaries take precedence over material interaction at a coincident hit.
Geometry is copied on insertion. Tolerance defaults to 1e-8 m; features thinner
than the tolerance are unresolved. `scene.validate()` validates primitives;
intersection relationships are resolved during tracing rather than by destructive
boolean edits to geometry.

The time convention is exp(-i omega t). A passive dielectric has
`epsilon_r = epsilon_real + 1j * epsilon_loss`. Conductivity in S/m contributes
an additional `+1j * sigma / (omega * EPS0)`. Do not count the same measured loss
in both epsilon's imaginary part and conductivity. Permeability is positive real.

## Gaussian-beam branching, fields and limits

Every dielectric interface deterministically produces reflected and transmitted
branches, except when transmission is non-propagating under total internal
reflection. Each branch stores its parent, interaction history, current region,
polarization, complex field, optical path, transverse frame and complex curvature.
Power transmission includes admittance and projected-area factors. Fields retain
Fresnel phase, propagation phase, absorption and continuous Gouy phase.

`result.field_at(position)` and `fields_at(positions)` sum the complex fields of
all segments covering the sample's longitudinal position. Field samples should be
away from interfaces; finite beam footprints at an interface are modeled locally.
Samples exactly at the final range/termination plane are excluded. Visualization
shows branch axes, not beam boundaries. `plot_field` plots sampled field magnitude.

`max_interactions`, `max_distance`, `max_branches` and `min_power_fraction`
control workload. Inspect `terminals` and `power_balance` for truncation; budget
termination is never relabeled as material absorption. Converge coherent results
with respect to depth, power cutoff and branch budget. A zero power cutoff disables
pruning of nonzero branches. `absorbed_power` includes propagation loss and opaque
metal absorption. `interface_flux_residual` is separate: in lossy incident media,
independent incident/reflected branch powers omit their interference flux. This
signed budget residual is not labeled as physical absorption.

## Beamlets and existing antennas

```python
beamlets = beam.beamlets(samples=25, radius_ratio=0.5)
result = rt.solve(scene, beamlets, config)
```

This is a coherent Gaussian-convolution decomposition; increase samples to check
spatial convergence. Narrower constituents resolve illumination spatially, but
must still satisfy the paraxial divergence limit. With multiple sources,
`power_balance` is the sum of constituent branch budgets, not physical coherent
power. Integrate the summed field flux to measure physical power.

`ApertureSource(position, direction, half_widths, samples, beamlet_radius, field)`
fits a callable transverse Cartesian complex field on a planar grid using Gaussian
basis functions. Its width, aperture extent and sample density require separate
convergence checks. `ApertureSource.from_antenna(...)` samples the existing native
`Antenna`/`load_farfield` models on a vacuum plane in their far zone; the solver
frequency must match that sampling frequency. Each aperture covers a forward
paraxial cone, not a full spherical antenna pattern.

## Gaussian-beam model scope

- Homogeneous isotropic, nondispersive regions; constant dielectric epsilon and
  permeability, with the conductivity term evaluated at the simulation frequency.
- Real Snell central-ray geometry with weak-loss attenuation. Default maximum
  Im(n)/Re(n) is 0.01; stronger dielectric loss is rejected. Normal-incidence
  lossy slab results are checked against the exact complex slab solution.
- Oblique lossy propagation is a weak-loss approximation. Lossy critical-angle
  events terminate with `unsupported_lossy_critical`; no fictitious transmitted
  ray is created. Strong-loss complex-ray transport and evanescent tunneling are
  not supported.
- Paraxial Gaussian propagation. Source divergence must be <=0.2 rad. Curved
  interfaces use local differential ray maps, including astigmatism. Wavelength
  and footprint must be small relative to relevant curvature/geometry scales.
- Whole-beam boundary interaction follows the central axis. Partial interception,
  clipping and edge diffraction are not solved exactly; footprint/edge proximity
  emits warnings. Beamlet refinement can resolve central-axis classification,
  but does not turn this into a full-wave diffraction solver.
- `save_result`/`load_result` use a new explicit JSON schema containing scene,
  configuration and complete branch states. Old HDF5/CSV results are not migrated.

## Inspect Gaussian-beam Tx/Rx patterns before tracing (VTK/PyVista)

```python
receiver_look = rt.GaussianBeam([0, 0, 1.0], [0, 0, -1], 0.04)
viewer = rt.inspect_scene(
    scene, transmitters=beam, receivers=receiver_look,
    frequency_hz=77e9,
)
viewer.close()
# Only now run rt.solve(...).
```

No tracing is performed by the inspector. Dielectrics are blue, opaque material
boundaries gold, termination boundaries wireframe, Tx red and Rx green. Gaussian
patterns show their 1/e field envelopes in the launch medium, plus their waist,
position and axis. Receiver direction means outward look direction. These Rx
descriptors are visualization inputs; they do not yet compute receive-port mode
overlap. `AntennaPattern(position, antenna)` displays normalized native/imported
antenna lobes with the antenna's orientation. Multiple Tx/Rx descriptors are
accepted as lists. `ApertureSource` can also preview constituent beamlets.

Use `ViewOptions(beam_length=..., pattern_scale=..., geometry_opacity=...)` to
adjust the display. `show=False, path="preview.png"` renders off screen and returns
a Plotter; close it after use. The same VTK geometry/pattern rendering is used by
`visualize(result, transmitters=..., receivers=...)` after tracing. Previews show
unscattered source/receive patterns, not beams propagated through the scene.
Termination boundaries are wireframe and excluded from automatic camera fitting
by default, so large enclosures do not hide the source geometry. Set
`ViewOptions(frame_termination=True)` to fit the complete enclosure, or
`show_termination=False` to hide it.

## Tests, benchmarks and examples

```powershell
python -m pytest -q
python benchmarks/beam_transport.py
python example1.py --no-show
python example2.py --no-show
```

Analytical CSV data, convergence plots, pass/fail tolerances and platform metadata
are saved under `benchmarks/results/beams`. See `benchmarks/BEAMS.md` for the
reference equations and what each convergence sweep establishes. The native
CTest suite and spreading benchmarks validate the point-ray numerical routines.

The Gaussian-beam validation exercises the shared C++ implementation with finite sources and plane-wave interface
limits. Run it to regenerate the local report and convergence figures; generated
results are ignored by Git. See [BEAMS.md](benchmarks/BEAMS.md) for independent
theory, tolerances and model scope.

The native spreading and circular-pipe caustic validations exercise the separate
point-ray solver. Their results do not validate Gaussian-beam transport. The main
room and cylinder examples use that point-ray API to produce impulse responses.
