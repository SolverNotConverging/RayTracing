# RayTracing 0.5

Coherent paraxial Gaussian beams, reflected/transmitted branch trees, and explicit
material volumes. Python owns scene topology and beam transport; C++ provides
validated geometric intersections and antenna far-field models.

This is a breaking API update. There are no compatibility aliases for the 0.4
point-receiver solver. Its native numerical routines remain internal regression
references, not the public propagation API.

## Installation

Requires CPython 3.14+, a C++20 compiler, CMake 4.2+, Ninja, Eigen3,
nlohmann_json and shared HDF5. NumPy and Matplotlib provide beam algebra and plots.
The internal native reference/persistence library still links HDF5.

On this development machine, the Git-ignored local shortcut initializes MSVC and
selects the Git-ignored CMakeUserPresets.json:

```powershell
.\install.ps1
```

Other Windows users should open an x64 Visual Studio Developer PowerShell,
activate their virtual environment, and configure their own vcpkg path:

```powershell
$env:CMAKE_GENERATOR = 'Ninja'
$env:CMAKE_ARGS = '-DCMAKE_TOOLCHAIN_FILE=C:/your/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows -DCMAKE_CXX_SCAN_FOR_MODULES=OFF'
python -m pip install .
```

Install `eigen3`, `nlohmann-json` and `hdf5` for `x64-windows` with vcpkg first.
For development, install `.[dev]`; `--no-build-isolation` can be used once the
build tools are installed in the virtual environment. On Linux/macOS, install
the corresponding native development libraries and configure CMAKE_PREFIX_PATH.

Keep `python/` unmarked as a **Sources Root** in PyCharm. Adding it to the run
configuration's `PYTHONPATH` makes the unbuilt source package shadow the installed
wheel, producing an import error for `_core`. If it is already marked, right-click
`python/` and choose **Mark Directory as → Unmark as Sources Root**, or disable
**Add source roots to PYTHONPATH** in the run configuration. Verify that
`python -c "import raytracing; print(raytracing.__file__)"` points into
`.venv`'s `site-packages/raytracing` directory.

## A focused beam through a lossy dielectric slab

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

## Geometry and materials

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

## Branching, fields and limits

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

## Model scope

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

## Inspect Tx/Rx patterns before tracing (VTK/PyVista)

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
CTest suite and older spreading benchmarks remain internal numerical references.

The example scripts and notebooks use the new API: example 1 is a lossy slab
with an enclosing termination boundary; example 2 is a focused reflected beam.
