# RayTracing

Python simulations backed by a C++20 numerical solver. The solver launches rays,
refines specular PEC paths to a point receiver, deduplicates them, calculates
Jacobian spreading, and reconstructs complex electric fields and coherent impulse
responses in a homogeneous, nondispersive medium. Matplotlib plots responses;
PyVista displays geometry, antenna patterns, paths and polarization.

## Build requirements

| Component | Requirement |
| --- | --- |
| Python | CPython 3.14 or newer; use 3.14 for the commands below |
| Compiler | C++20-capable MSVC on Windows, Apple Clang on macOS, or GCC on Linux |
| Build tools | CMake 4.2+, Ninja, scikit-build-core and pybind11 3.0+ |
| Native libraries | Eigen3, nlohmann_json and **shared** HDF5 |
| Python runtime | NumPy, Matplotlib and PyVista; installed with the package |
| Notebook tools | JupyterLab and ipykernel; installed with the `notebooks` extra |

The numerical core links Eigen3. The persistence layer additionally uses
nlohmann_json to serialize simulation records inside HDF5 and links shared HDF5.
Both are required for the complete Python package. Rendering is implemented in
Python; PyVista installs its own VTK dependency.

Windows/MSVC builds, both examples and the regression suites have been exercised
on the development machine. The macOS and Linux recipes below follow the current
CMake configuration and dependency requirements; they have not yet been tested
on those operating systems.

## Prepare the source

Extract a source release or obtain a repository checkout, then open a terminal
in its root directory, next to `pyproject.toml` and `CMakeLists.txt`.

Use the existing `.venv` if present. The creation commands below are for a fresh
checkout; skip them when reusing an environment. Keep Python, compiler and native
libraries on the same CPU architecture.

## Windows: MSVC

Install 64-bit CPython 3.14, Git, and Visual Studio or Build Tools with the
**Desktop development with C++** workload, including MSVC and a Windows SDK.
Open **Developer PowerShell for Visual Studio** targeting x64, then change into
the source root. `cl` must be available in this shell.

```powershell
# Create only if .venv does not already exist.
py -3.14 -m venv .venv
. ./.venv/Scripts/Activate.ps1
python -m ensurepip --upgrade
python -m pip install --upgrade pip
python -m pip install "scikit-build-core>=0.11" "pybind11>=3.0" "cmake>=4.2" ninja build

# Use your existing vcpkg checkout, or bootstrap one at this location.
$env:VCPKG_ROOT = "$env:USERPROFILE/vcpkg"
git clone https://github.com/microsoft/vcpkg.git $env:VCPKG_ROOT
& "$env:VCPKG_ROOT/bootstrap-vcpkg.bat"
$env:VCPKG_DEFAULT_TRIPLET = 'x64-windows'
& "$env:VCPKG_ROOT/vcpkg.exe" install eigen3 nlohmann-json hdf5 --triplet $env:VCPKG_DEFAULT_TRIPLET

$env:CMAKE_GENERATOR = 'Ninja'
$env:CMAKE_ARGS = "-DCMAKE_TOOLCHAIN_FILE=`"$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake`" -DVCPKG_TARGET_TRIPLET=$env:VCPKG_DEFAULT_TRIPLET -DCMAKE_CXX_SCAN_FOR_MODULES=OFF"
python -m pip install --no-build-isolation ".[dev,notebooks]"
```

If vcpkg already exists, set `VCPKG_ROOT` to that directory and skip the clone.
If PowerShell blocks activation, invoke `.venv/Scripts/python.exe` explicitly
and add `.venv/Scripts` to this shell's `PATH` for CMake and Ninja. Keep native
builds in the developer shell so the compiler, linker, headers and SDK are found.
See Microsoft's [vcpkg host prerequisites](https://learn.microsoft.com/vcpkg/concepts/supported-hosts).

## macOS: Apple Clang and Homebrew

Install Xcode Command Line Tools and [Homebrew](https://brew.sh/). Use a native
arm64 terminal/Python on Apple Silicon, or x86_64 on Intel; keep the architectures
consistent. Install Python and the C++ libraries with Homebrew:

```bash
xcode-select --install  # Skip if Command Line Tools are already installed.
brew install python@3.14 eigen@3 nlohmann-json hdf5

# From the RayTracing source root; create only if .venv does not exist.
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install 'scikit-build-core>=0.11' 'pybind11>=3.0' 'cmake>=4.2' ninja build

export CC="$(xcrun --find clang)"
export CXX="$(xcrun --find clang++)"
export CMAKE_PREFIX_PATH="$(brew --prefix eigen@3):$(brew --prefix nlohmann-json):$(brew --prefix hdf5)"
export CMAKE_GENERATOR=Ninja
export CMAKE_ARGS="-DCMAKE_CXX_SCAN_FOR_MODULES=OFF -DCMAKE_INSTALL_RPATH=\"$(brew --prefix hdf5)/lib\""
python -m pip install --no-build-isolation '.[dev,notebooks]'
```

The explicit prefixes let CMake find both Apple Silicon and Intel Homebrew
installations, including the versioned [Eigen 3 formula](https://formulae.brew.sh/formula/eigen@3).
The extension links Homebrew's [HDF5](https://formulae.brew.sh/formula/hdf5);
keep that package installed when using this local build.

## Linux: GCC and distribution packages

Use a C++20-capable GCC toolchain (GCC 12+ is a suitable starting point). Install
the development headers and libraries through your distribution's package manager.
For Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install build-essential pkg-config libeigen3-dev nlohmann-json3-dev libhdf5-dev
```

For Fedora:

```bash
sudo dnf install gcc gcc-c++ make pkgconf-pkg-config eigen3-devel json-devel hdf5-devel
```

Install CPython 3.14 using your distribution or a Python installer. If using the
distribution's Python, also install its matching development and venv packages
(for example `python3.14-dev`/`python3.14-venv` on Debian/Ubuntu when available,
or the matching `python3-devel` on Fedora).
If your distribution does not provide Python 3.14,
[uv](https://docs.astral.sh/uv/guides/install-python/) can supply it with
`uv python install 3.14` and create the environment with
`uv venv --python 3.14 .venv`. Otherwise:

```bash
# From the RayTracing source root; skip creation if .venv already exists.
python3.14 -m venv .venv
```

Then activate that environment and build against the system libraries:

```bash
source .venv/bin/activate
python -m ensurepip --upgrade
python -m pip install --upgrade pip
python -m pip install 'scikit-build-core>=0.11' 'pybind11>=3.0' 'cmake>=4.2' ninja build

export CC=gcc
export CXX=g++
export CMAKE_GENERATOR=Ninja
export CMAKE_ARGS='-DCMAKE_CXX_SCAN_FOR_MODULES=OFF'
python -m pip install --no-build-isolation '.[dev,notebooks]'
```

CMake discovers Eigen3 and nlohmann_json through their installed package configs.
Its standard `FindHDF5` module finds the system HDF5 headers and libraries,
including distribution-specific serial-library layouts. Keep the HDF5 runtime
package installed when using this local build.

All three recipes disable C++ module scanning: this project uses headers rather
than C++ modules, so a separate Clang dependency scanner is unnecessary.
See CMake's [module-scanning setting](https://cmake.org/cmake/help/latest/variable/CMAKE_CXX_SCAN_FOR_MODULES.html).

## Verify installation and run examples

With `.venv` active, these commands are the same on all three operating systems:

```bash
python -c "import raytracing as rt; print(rt.__version__); print(rt.example_pattern())"
python example1.py
python example2.py
```

Each script solves the scene, saves HDF5/CSV/PNG files under `results/exampleN`,
and opens an interactive desktop window. Left-drag rotates, Shift+drag or
middle-drag pans, the wheel zooms, and **R** resets the camera. Close the first
window to finish its script before starting the second.

| Example | Scene |
| --- | --- |
| [example1.py](example1.py) / [example1.ipynb](example1.ipynb) | Mixed geometry, CST transmitter aimed horizontally along +x, isotropic receiver |
| [example2.py](example2.py) / [example2.ipynb](example2.ipynb) | Cylinder of radius 0.127 m from z=0 to -0.5 m, open top and bottom disk; both CST patches face -z at (0.006, 0, 0) m |

Example 2's common 6 mm lateral shift avoids the exact axial caustic while
retaining monostatic reception. Set `source_shift = 0` to inspect that singular case.

For notebooks:

```bash
python -m ipykernel install --sys-prefix --name raytracing --display-name "RayTracing (.venv)"
python -m jupyterlab
```

Select the **RayTracing (.venv)** kernel. Notebooks display plots and scene images
inline. Use the `.py` examples for desktop interaction. Restart an existing kernel
after reinstalling the package.

In PyCharm, select the existing interpreter at `.venv/Scripts/python.exe` on
Windows or `.venv/bin/python` on macOS/Linux, then run the example as a Python
script. `tool.uv.managed = false` keeps `uv run` from replacing the installed
wheel with an automatic editable build. Explicit installation commands manage
the environment. Changing a simulation script needs no rebuild; changing the
library sources requires a reinstall.

Keep `python/` unmarked as a **Sources Root** in PyCharm. Adding it to the run
configuration's `PYTHONPATH` makes the unbuilt source package shadow the installed
wheel, producing an import error for `_core`. If it is already marked, right-click
`python/` and choose **Mark Directory as → Unmark as Sources Root**, or disable
**Add source roots to PYTHONPATH** in the run configuration. Verify that
`python -c "import raytracing; print(raytracing.__file__)"` points into
`.venv`'s `site-packages/raytracing` directory.

Desktop viewing requires a graphical session with working OpenGL support. For
remote/headless simulation, call `rt.solve` and save results without opening a
viewer. `rt.plot` can use Matplotlib's `Agg` backend. PyVista screenshots still
require a functioning rendering backend even with `off_screen=True`.

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
| `show(result, options=None)` | Open an interactive desktop viewer |

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

Open an interactive desktop window directly from a standalone Python script:

```python
rt.show(saved)  # Or rt.show(result) immediately after solving.
```

Drag with the left mouse button to rotate, Shift+drag or middle-drag to pan,
and scroll to zoom. Press **R** to reset the camera. Close the window to resume
the script. `show` accepts the same `ViewOptions` for path
selection and appearance, and automatically enables desktop interaction.

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

## Rebuild, test and share

Keep the compiler environment and `CMAKE_ARGS` from your platform's setup active.
After editing the package sources, close Python processes that have loaded the
native extension, then rebuild and reinstall in the same `.venv`:

```bash
python -m pip install --no-build-isolation --no-deps --force-reinstall .
```

The normal install builds a Release extension. `--no-deps` here preserves the
already installed dependencies; omit it when dependency requirements change.
If you change the compiler, architecture or Python version, use a fresh CMake
build directory rather than reusing a cache from the previous toolchain.

### Tests

The `dev` extra installed above includes pytest. Run the binding, persistence
and rendering checks with:

```bash
python -m pytest tests/python
```

Rendering tests need a working graphics backend. For the seven C++ suites,
configure a separate build in the source root. Windows Developer PowerShell:

```powershell
cmake -S . -B cmake-build-tests -G Ninja -DCMAKE_BUILD_TYPE=Release `
  "-DCMAKE_TOOLCHAIN_FILE=$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" `
  "-DVCPKG_TARGET_TRIPLET=$env:VCPKG_DEFAULT_TRIPLET" `
  "-DPython_EXECUTABLE=$((Get-Command python).Source)" `
  -DCMAKE_CXX_SCAN_FOR_MODULES=OFF -DBUILD_TESTING=ON -DBUILD_BENCHMARKS=ON
cmake --build cmake-build-tests --parallel
ctest --test-dir cmake-build-tests --output-on-failure
```

macOS/Linux, using the compiler and environment from the package build
(including `CMAKE_PREFIX_PATH` on macOS):

```bash
cmake -S . -B cmake-build-tests -G Ninja -DCMAKE_BUILD_TYPE=Release \
  "-DPython_EXECUTABLE=$(command -v python)" \
  -DCMAKE_CXX_SCAN_FOR_MODULES=OFF -DBUILD_TESTING=ON -DBUILD_BENCHMARKS=ON
cmake --build cmake-build-tests --parallel
ctest --test-dir cmake-build-tests --output-on-failure
```

CMake exposes `RayTracing::RayTracing` for numerics and `RayTracing::IO` for
persistence. Python bindings live in `python/bindings`; rendering lives in
`python/raytracing/visualization.py`.

The suites cover surfaces, sequence replay, refinement, fields, antennas and
persistence, cylinder returns and analytical spreading. Python checks cover
NumPy ownership, pattern rotation, shifted cylinder returns, saved results and
rendering. See [benchmarks/README.md](benchmarks/README.md) for the analytical
formulas and Matplotlib comparison plots.

### Distributions

Build a wheel and a source archive without creating another Python environment:

```bash
python -m build --wheel --sdist --no-isolation
```

Outputs appear in `dist/`. The source archive includes the CST patch, Python/C++
sources and both script/notebook examples. A recipient extracts it and follows
the matching platform instructions above.

A wheel is specific to its operating system, CPU architecture and Python ABI.
For example, `raytracing-0.3.0-cp314-cp314-win_amd64.whl` requires CPython 3.14 on
64-bit Windows. Recipients install a matching wheel in their active environment:

```bash
python -m pip install /path/to/raytracing-0.3.0-cp314-cp314-win_amd64.whl
```

Replace that example filename with the wheel built for the recipient. Installing
a wheel does not require a compiler. Share the example scripts/notebooks alongside
it or provide the source archive; they are not installed as top-level scripts by
the wheel. The wheel includes the CST file accessed through `rt.example_pattern()`.

Windows wheels produced by this project bundle HDF5 and zlib DLLs. The recipient
also needs the Microsoft Visual C++ runtime for their architecture. macOS/Linux
wheels from the recipes above use Homebrew or system shared libraries. Before distributing them as standalone binaries, bundle their native
dependencies with [delocate](https://github.com/matthew-brett/delocate) on macOS or
[auditwheel](https://github.com/pypa/auditwheel) on Linux, and test the repaired
wheel on a clean target machine. For broad Linux compatibility, build in an
appropriate manylinux environment; repairing a wheel does not lower the compiler
or glibc requirements of code already built.

### Troubleshooting

| Symptom | Check or fix |
| --- | --- |
| CMake cannot find a C/C++ compiler | Windows: use Developer PowerShell with the C++ workload installed. macOS: install Command Line Tools. Linux: install GCC/G++. |
| CMake cannot find Eigen3, nlohmann_json or HDF5 | Windows: check the vcpkg toolchain/triplet. macOS: check the Homebrew packages and `CMAKE_PREFIX_PATH`. Linux: install the development packages listed above. |
| Import fails with an HDF5 `.so` or `.dylib` error | Restore the Homebrew/system HDF5 runtime and rebuild if its location or ABI changed, or use a repaired wheel with bundled dependencies. |
| Windows import fails with a DLL error | Install the Visual C++ runtime and confirm that `raytracing/.libs` contains the wheel's HDF5/zlib DLLs. |
| Build reports a missing `simple_patch.ffs` | Restore `antenna_patterns/simple_patch.ffs` from the source release before installing. |
| A new library function is missing | Confirm the selected `.venv`, reinstall the package and restart the notebook kernel or Python process. |
| IDE tries to build before running a script | Use the existing `.venv` interpreter directly, or `uv run --no-sync`; keep `tool.uv.managed = false`. |
| Notebook desktop window does not respond | Run `example1.py` or `example2.py` as a standalone script; use inline previews in notebooks. |
| Python loads standard-library names from a cache directory | Remove that directory from `PYTHONPATH` and exclude dependency caches/build folders from IDE source roots. |
