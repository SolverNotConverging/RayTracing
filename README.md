# RayTracing

A C++20 electromagnetic ray-tracing library with a Python API. It contains two solvers:

- **Point-ray multipath** (`rt::solve`): launches rays from a transmitter, refines exact
  specular paths to a point receiver, computes geometrical spreading (with caustic phase)
  and Fresnel reflections, and returns a coherent impulse response.
- **Gaussian beams** (`rt::beam::solve`): paraxial beam transport through dielectric
  volumes and conductor/PEC surfaces, with reflection, refraction and loss.

Both share the same geometry, materials and antenna models. SI units throughout;
phasors use the `exp(-i omega t)` convention.

## Python

Requires Python 3.14+, a C++20 compiler, CMake 4.2+ and Eigen3.

```bash
pip install .[viz]
```

On Windows with vcpkg, run from an x64 Developer PowerShell and pass the toolchain, for example
`$env:CMAKE_ARGS = '-DCMAKE_TOOLCHAIN_FILE=C:/opt/vcpkg/scripts/buildsystems/vcpkg.cmake'`.

```python
import raytracing as rt

scene = rt.Scene()
scene.add(rt.Rectangle([2, 0, 0], [0, 1, 0], [0, 0, 1], 3, 1), rt.Material.pec())

tx = rt.Transmitter([0, 0, 0], rt.ThinWireDipole(rt.SPEED_OF_LIGHT / 77e9 / 2))
rx = rt.Receiver([1, 0.5, 0], rt.Isotropic())
result = rt.solve(scene, tx, rx, rt.SolverConfig(frequency_hz=77e9))

print(result.impulse_response.delays, result.impulse_response.coefficients)
rt.show_paths(result)
```

```python
from raytracing import beam

scene = beam.Scene()
scene.add_volume(rt.Box([0, 0, 0.5], [1, 1, 0.01]), rt.Material.dielectric(4 + 0.02j))
result = beam.solve(scene, beam.GaussianBeam(waist_radius=0.1), beam.SolverConfig(frequency_hz=10e9))
print(result.field_at([0, 0, 1]))
rt.show_beams(result)
```

`plot_impulse_response`, `plot_field`, `show_paths`, `show_beams` and `preview` need the
`viz` extra (Matplotlib, PyVista). Results expose NumPy arrays, so save what you need with
`numpy.savez`.

## C++

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
cmake --install build --prefix <prefix>
```

Consume it with `find_package` (set `CMAKE_PREFIX_PATH` to the install prefix) or
`add_subdirectory`:

```cmake
find_package(RayTracing 1 CONFIG REQUIRED)
target_link_libraries(app PRIVATE RayTracing::RayTracing)
```

```cpp
#include <raytracing/solver.hpp>

rt::Scene scene;
scene.add(rt::Rectangle{{2, 0, 0}, {0, 1, 0}, {0, 0, 1}, 3, 1}, rt::Material::pec());
rt::Transmitter tx{{0, 0, 0}, rt::ThinWireDipole{rt::speedOfLight / 77e9 / 2}};
rt::Receiver rx{{1, 0.5, 0}, rt::Isotropic{}};
const auto result = rt::solve(scene, tx, rx);
for (const auto &tap: result.impulseResponse.taps) std::cout << tap.delay << ' ' << tap.coefficient << '\n';
```

`<raytracing/raytracing.hpp>` includes everything. See [examples](examples) for complete
C++ and Python programs.

## Model notes

- **Materials** are `dielectric(epsilon, conductivity, mu)`, `conductor(conductivity)` or `pec()`.
  Point-ray surfaces reflect as the boundary of a half-space of their material; beam
  dielectrics must be volumes (`Box`, `Sphere` or capped `Cylinder`).
- **Antennas**: `Isotropic`, `ShortDipole`, `ThinWireDipole`, `RectangularAperture`, or a
  CST `.ffs` / HFSS `.ffd` far field via `load_farfield`. Receive responses are reciprocal.
- **Validation** happens once, when shapes and materials are added and when a solver starts.
  The low-level routines (`intersect`, `refine_path`, ...) assume valid input.
- **Beam limits**: paraxial (divergence <= 0.2 rad), weak dielectric loss
  (`max_loss_ratio`), whole-beam interactions along the axis; edge clipping and
  diffraction are not modelled, only warned about.
- **Point-ray paths** whose spreading is singular (receiver caustics) keep their geometry
  but have no field.

## Layout

```
include/raytracing/   public headers          src/            library sources
python/bindings/      pybind11 module         python/raytracing/  Python package
tests/cpp, tests/python                       examples/cpp, examples/python
```
