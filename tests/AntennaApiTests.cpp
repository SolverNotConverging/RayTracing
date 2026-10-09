#include <RayTracing/ResultIO.hpp>
#include <RayTracing/Solver.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool v, const char *m) {
  if (!v)
    throw std::runtime_error(m);
}
void near(Complex a, Complex b, const char *m, double tolerance = 1e-9) {
  check(std::abs(a - b) < tolerance, m);
}
template <class F> void rejects(F f) {
  bool failed = false;
  try {
    f();
  } catch (const std::exception &) {
    failed = true;
  }
  check(failed, "Invalid antenna or data accepted");
}
void write_ffd(const std::filesystem::path &p, bool dependent) {
  std::ofstream out(p);
  out << "0 180 3\n0 360 5\n";
  if (dependent)
    out << "Frequencies 2\n";
  for (int f = 0; f < (dependent ? 2 : 1); ++f) {
    if (dependent)
      out << "Frequency " << (f == 0 ? 76e9 : 78e9) << '\n';
    for (int th = 0; th < 3; ++th)
      for (int ph = 0; ph < 5; ++ph)
        out << (th == 1 ? -(1 + 2 * f) : 0) << ' '
            << (th == 1 ? -(2 + 2 * f) : 0) << " 0 0\n";
  }
}
} // namespace
int main() try {
  using namespace rt;
  const auto folder = std::filesystem::path("antenna-api-test-results");
  std::filesystem::create_directories(folder);
  Antenna iso = Isotropic{};
  iso.validate(77e9);
  near(iso.farfield(Vec3::UnitX(), 77e9).norm(), 1,
       "Isotropic amplitude incorrect");
  Antenna small = ShortDipole{1e-4, 1};
  small.validate(77e9);
  const double expected = 376.730313668 * 2 * PI * 77e9 / C * 1e-4 / (4 * PI);
  near(small.farfield(Vec3::UnitX(), 77e9)[2], Complex(0, expected),
       "Short dipole analytical field incorrect");
  near(small.farfield(Vec3::UnitZ(), 77e9).norm(), 0, "Dipole axial null lost");
  near(std::abs(
           small.receive(Vec3::UnitX(), Vec3::UnitZ().cast<Complex>(), 77e9)),
       std::sqrt(1.5), "Short dipole receive gain incorrect");
  Antenna half = ThinWireDipole{C / 77e9 / 2, 1};
  half.validate(77e9);
  near(half.farfield(Vec3::UnitX(), 77e9)[2],
       Complex(0, 376.730313668 / (2 * PI)),
       "Half-wave dipole broadside field incorrect");
  Antenna aperture =
      RectangularAperture{2 * C / 77e9, 2 * C / 77e9, 1, 0, true};
  aperture.validate(77e9);
  check(aperture.farfield(Vec3::UnitZ(), 77e9).norm() > 0,
        "Aperture broadside missing");
  near(aperture.farfield(-Vec3::UnitZ(), 77e9).norm(), 0,
       "PEC-backed aperture radiates backwards");
  near(aperture.farfield(Vec3(0.5, 0, std::sqrt(0.75)), 77e9).norm(), 0,
       "Rectangular aperture sinc null missing");
  small.orientation =
      Eigen::AngleAxisd(PI / 2, Vec3::UnitY()).toRotationMatrix();
  near(small.farfield(Vec3::UnitZ(), 77e9)[0], Complex(0, expected),
       "Antenna orientation not applied");
  rejects([&] {
    auto bad = iso;
    bad.orientation(0, 0) = 2;
    bad.validate(77e9);
  });

  const auto patch = load_farfield(TEST_PATTERN);
  check(patch.pattern->frequenciesHz.size() == 5 &&
            patch.pattern->thetaDegrees.size() == 181 &&
            patch.pattern->phiDegrees.size() == 361,
        "CST grid dimensions incorrect");
  const auto value = patch.farfield(Vec3::UnitX(), 77e9);
  near(value[2], Complex(1.27639257, 3.28887389),
       "CST Etheta, phase convention or frequency ordering incorrect", 1e-7);
  near(value[1], Complex(-1.46396372e-4, 5.07387899e-4),
       "CST Ephi not preserved", 1e-9);
  near(patch.reference_power(77e9), 0.4978351,
       "CST accepted power metadata incorrect", 1e-8);
  auto rotatedPatch = patch;
  rotatedPatch.rotate(Vec3(0, 2, 0), -90.0); // Axis need not be normalized.
  const auto rotation = Eigen::AngleAxisd(-PI / 2, Vec3::UnitY()).toRotationMatrix();
  check((rotatedPatch.orientation * Vec3::UnitZ() + Vec3::UnitX()).norm() < 1e-12,
        "Degree rotation did not turn local z horizontal");
  const Vec3 probe = Vec3(1, 2, 3).normalized();
  check((rotatedPatch.farfield(rotation * probe, 77e9) -
         rotation.cast<Complex>() * patch.farfield(probe, 77e9)).norm() < 1e-10,
        "Rotation did not rotate imported field direction and polarization together");
  const Vec3C incident = Vec3(2, -1, 0).cast<Complex>();
  near(rotatedPatch.receive(rotation * probe, rotation.cast<Complex>() * incident, 77e9),
       patch.receive(probe, incident, 77e9), "Rx rotation changed reciprocal contraction");
  rotatedPatch.rotate(Vec3::UnitZ(), 30.0);
  const Eigen::Matrix3d composed = Eigen::AngleAxisd(PI / 6, Vec3::UnitZ()).toRotationMatrix() * rotation;
  check((rotatedPatch.orientation - composed).norm() < 1e-12,
        "Rotation did not compose about world axis");
  rotatedPatch.validate(77e9);
  rejects([&] { rotatedPatch.rotate(Vec3::Zero(), 90.0); });
  rejects([&] { rotatedPatch.rotate(Vec3::UnitY(), std::numeric_limits<double>::infinity()); });
  check((rotatedPatch.orientation - composed).norm() < 1e-12,
        "Invalid rotation changed antenna orientation");
  const auto mid = patch.farfield(Vec3::UnitX(), 77.25e9);
  const auto high = patch.farfield(Vec3::UnitX(), 77.5e9);
  check((mid - (value + high) / 2).norm() < 1e-10,
        "Complex frequency interpolation incorrect");
  rejects([&] { patch.farfield(Vec3::UnitX(), 80e9); });
  const auto nearPole = patch.farfield(Vec3(1e-9, 0, 1).normalized(), 77e9);
  const auto otherPole = patch.farfield(Vec3(0, 1e-9, 1).normalized(), 77e9);
  check((nearPole - otherPole).norm() < 1e-6,
        "Cartesian interpolation discontinuous at pole");

  write_ffd(folder / "dependent.ffd", true);
  write_ffd(folder / "independent.ffd", false);
  FarfieldImportOptions import;
  import.inputConvention = PhasorConvention::NegativeTime;
  import.inputPowerWatts = 1;
  const auto hfss = load_farfield(folder / "dependent.ffd", import);
  near(hfss.farfield(Vec3::UnitX(), 77e9)[2], Complex(2, 3),
       "HFSS block ordering or complex interpolation incorrect");
  const auto independent = load_farfield(folder / "independent.ffd", import);
  near(independent.farfield(Vec3::UnitX(), 50e9)[2], Complex(1, 2),
       "Frequency-independent HFSS file incorrect");
  const auto noPower = load_farfield(folder / "dependent.ffd");
  rejects([&] {
    noPower.receive(Vec3::UnitX(), Vec3::UnitZ().cast<Complex>(), 77e9);
  });
  // Reciprocal receive contraction is bilinear for a transmit-pattern
  // representation.
  near(hfss.receive(Vec3::UnitX(),
                    Complex(0, 1) * Vec3::UnitZ().cast<Complex>(), 77e9),
       Complex(-3, 2) * std::sqrt(4 * PI / (2 * 376.730313668)),
       "Receive pattern accidentally conjugated");

  SolverConfig cfg;
  cfg.rayCount = 64;
  cfg.maxReflections = 0;
  cfg.maxDistance = 5;
  Scene scene;
  Transmitter tx{{0, 0, 0}, Isotropic{}};
  Receiver rx{{3, 0, 0}, Isotropic{}};
  const auto direct = solve(scene, tx, rx, cfg);
  check(direct.rays.size() == 1 && direct.rays[0].field,
        "Solver did not return unique refined direct ray");
  near(std::abs(direct.impulseResponse.taps_[0].coefficient_), 1.0 / 3,
       "API direct field gain incorrect");
  tx.antenna = Isotropic{Polarization::Vertical, 2};
  const auto scaled = solve(scene, tx, rx, cfg);
  near(scaled.impulseResponse.taps_[0].coefficient_,
       2.0 * direct.impulseResponse.taps_[0].coefficient_,
       "Source directional amplitude normalized away");
  tx.antenna = ShortDipole{};
  rx.position = Vec3(0, 0, 3);
  const auto null = solve(scene, tx, rx, cfg);
  check(null.rays.size() == 1 && null.rays[0].field,
        "Source null prevented geometry discovery");
  near(null.impulseResponse.taps_[0].coefficient_, 0,
       "Source null generated field");

  // Colocated patch antennas must find a physical return rather than a direct ray.
  Scene tube;
  tube.add(Cylinder{{0, 0, -0.25}, Vec3::UnitZ(), 0.127, 0.25, false});
  tube.add(Disk{{0, 0, -0.5}, -Vec3::UnitZ(), 0.127});
  auto monostaticConfig = cfg;
  monostaticConfig.rayCount = 512;
  monostaticConfig.maxReflections = 8;
  const Transmitter monostaticTx{Vec3::Zero(), patch};
  const Receiver monostaticRx{Vec3::Zero(), patch};
  const auto returned = solve(tube, monostaticTx, monostaticRx, monostaticConfig);
  check(returned.launchedRays == monostaticConfig.rayCount,
        "Colocated solver launched a zero-length direct direction");
  check(!returned.impulseResponse.taps_.empty(), "Colocated patch return missing");
  bool bottomReturn = false;
  for (const auto &ray : returned.rays) {
    check(ray.refinement.geometry_.receiver_->pathDistance_ > 0,
          "Colocated solver returned a zero-length path");
    const auto &hits = ray.refinement.geometry_.reflections_;
    if (hits.size() == 1 && hits[0].surfaceIndex_ == 1) {
      bottomReturn = true;
      near(ray.refinement.geometry_.receiver_->pathDistance_, 1.0,
           "Bottom return should travel one metre", 1e-6);
      check(ray.field && std::isfinite(std::abs(ray.receivedCoefficient)) &&
                std::abs(ray.receivedCoefficient) > 0,
            "Colocated patch return has no finite field");
    }
  }
  check(bottomReturn, "Bottom disk reflection was not refined");
  const auto noReturn = solve(Scene{}, monostaticTx, monostaticRx, monostaticConfig);
  check(noReturn.rays.empty() && noReturn.impulseResponse.taps_.empty(),
        "Empty monostatic scene generated a direct return");

  // Embed the real CST pattern at both ends, plus every geometry variant.
  Scene full;
  full.add(Sphere{{10, 0, 0}, 1});
  full.add(Disk{{11, 0, 0}, Vec3::UnitX(), 1});
  full.add(Rectangle{{12, 0, 0}, Vec3::UnitY(), Vec3::UnitZ(), 1, 1});
  full.add(Cylinder{{13, 0, 0}, Vec3::UnitZ(), 1, 1, true});
  full.add(Triangle{{14, 0, 0}, {14, 1, 0}, {14, 0, 1}});
  tx.antenna = patch;
  rx.position = Vec3(3, 0, 0);
  rx.antenna = rotatedPatch;
  auto patterned = solve(full, tx, rx, cfg);
  save_h5(patterned, folder / "simulation.h5");
  const auto loaded = load_h5(folder / "simulation.h5");
  check(loaded.scene.surfaces().size() == 5 &&
            loaded.rays.size() == patterned.rays.size() &&
            loaded.candidates.size() == patterned.candidates.size(),
        "HDF5 geometry or path state lost");
  check(loaded.settings.rayCount == cfg.rayCount &&
            loaded.launchedRays == patterned.launchedRays,
        "HDF5 settings/diagnostics lost");
  check((loaded.receiver.antenna.orientation - rotatedPatch.orientation).norm() < 1e-12,
        "HDF5 rotated antenna orientation lost");
  near(loaded.impulseResponse.taps_[0].coefficient_,
       patterned.impulseResponse.taps_[0].coefficient_,
       "HDF5 response changed");
  check((loaded.rays[0].field->receiverField_ -
         patterned.rays[0].field->receiverField_)
                .norm() < 1e-12,
        "HDF5 incident field changed");
  check((loaded.transmitter.antenna.farfield(Vec3::UnitX(), 77e9) - value)
                .norm() < 1e-12,
        "HDF5 embedded pattern changed");
  auto changed = loaded;
  update_receiver(changed, Isotropic{});
  check(changed.rays[0].field->receiverField_ ==
            loaded.rays[0].field->receiverField_,
        "Changing Rx retraced/changed incident field");
  check(std::abs(changed.rays[0].receivedCoefficient -
                 loaded.rays[0].receivedCoefficient) > 1e-4,
        "Rx pattern change had no effect");
  save_impulse_csv(patterned.impulseResponse, folder / "response.csv");
  const auto csv = load_impulse_csv(folder / "response.csv");
  near(csv.taps_[0].coefficient_,
       patterned.impulseResponse.taps_[0].coefficient_,
       "CSV complex coefficient lost");
  check(csv.taps_[0].pathIndices_ ==
                patterned.impulseResponse.taps_[0].pathIndices_ &&
            csv.frequencyHz_ == 77e9,
        "CSV metadata lost");
  save_h5(null, folder / "null.h5");
  near(load_h5(folder / "null.h5").impulseResponse.taps_[0].coefficient_, 0,
       "Zero field HDF5 round trip failed");
  cfg.maxDistance = 0.1;
  const auto empty = solve(scene, tx, rx, cfg);
  save_h5(empty, folder / "empty.h5");
  check(load_h5(folder / "empty.h5").rays.empty(),
        "Empty result HDF5 round trip failed");
  save_impulse_csv(empty.impulseResponse, folder / "empty.csv");
  check(load_impulse_csv(folder / "empty.csv").taps_.empty(),
        "Empty response CSV failed");
  std::cout << "Antenna models, CST/HFSS imports, solver API and HDF5/CSV "
               "round trips passed\n";
  return 0;
} catch (const std::exception &e) {
  std::cerr << e.what() << '\n';
  return 1;
}
