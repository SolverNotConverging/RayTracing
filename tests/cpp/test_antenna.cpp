#include <raytracing/solver.hpp>
#include <raytracing/core.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace rt;

namespace {
constexpr double C = speedOfLight;
constexpr double PI = pi;
constexpr double ETA0 = freeSpaceImpedance;

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
  const auto folder = std::filesystem::path("antenna-api-test-results");
  std::filesystem::create_directories(folder);
  Antenna iso = Isotropic{};
  iso.validate(77e9);
  near(iso.farfield(Vec3::UnitX(), 77e9).norm(), 1,
       "Isotropic amplitude incorrect");
  Antenna small = ShortDipole{1e-4, 1};
  small.validate(77e9);
  const double expected = ETA0 * 2 * PI * 77e9 / C * 1e-4 / (4 * PI);
  near(small.farfield(Vec3::UnitX(), 77e9)[2], Complex(0, expected),
       "Short dipole analytical field incorrect");
  near(small.farfield(Vec3::UnitZ(), 77e9).norm(), 0, "Dipole axial null lost");
  near(std::abs(
           small.receive(Vec3::UnitX(), Vec3::UnitZ().cast<Complex>(), 77e9)),
       std::sqrt(1.5), "Short dipole receive gain incorrect");
  Antenna half = ThinWireDipole{C / 77e9 / 2, 1};
  half.validate(77e9);
  near(half.farfield(Vec3::UnitX(), 77e9)[2],
       Complex(0, ETA0 / (2 * PI)),
       "Half-wave dipole broadside field incorrect");
  {
    // R_rad = eta Cin(2 pi) / (4 pi); Cin(x) = int_0^x (1 - cos t) / t dt by
    // composite Simpson, independent of the antenna quadrature.
    const int n = 200000;
    const double h = 2 * PI / n;
    double cin = 0;
    for (int i = 0; i <= n; ++i) {
      const double t = i * h;
      const double f = i == 0 ? 0.0 : (1 - std::cos(t)) / t;
      cin += (i == 0 || i == n ? 1 : (i % 2 ? 4 : 2)) * f;
    }
    cin *= h / 3;
    const double expectedPower = 0.5 * ETA0 * cin / (4 * PI); // |I| = 1 A
    check(std::abs(half.reference_power(77e9) / expectedPower - 1) < 1e-9,
          "Half-wave dipole radiated power quadrature inaccurate");
  }
  Antenna aperture =
      RectangularAperture{2 * C / 77e9, 2 * C / 77e9, 1, 0};
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

  const auto patch = load_farfield(TEST_DATA "/simple_patch.ffs");
  const auto &grid = *std::get<std::shared_ptr<const FarfieldPattern>>(patch.model);
  check(grid.frequenciesHz.size() == 5 && grid.thetaDegrees.size() == 19 &&
            grid.phiDegrees.size() == 37,
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
  import.convention = PhasorConvention::NegativeTime;
  import.inputPower = 1;
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
       Complex(-3, 2) * std::sqrt(4 * PI / (2 * ETA0)),
       "Receive pattern accidentally conjugated");

  SolverConfig cfg;
  cfg.rayCount = 64;
  cfg.maxReflections = 0;
  cfg.maxDistance = 5;
  Scene scene;
  Transmitter tx{{0, 0, 0}, Isotropic{}};
  Receiver rx{{3, 0, 0}, Isotropic{}};
  const auto direct = solve(scene, tx, rx, cfg);
  check(direct.paths.size() == 1 && direct.paths[0].field,
        "Solver did not return unique refined direct ray");
  near(std::abs(direct.impulseResponse.taps[0].coefficient), 1.0 / 3,
       "API direct field gain incorrect");
  check(direct.candidates.size() == 1, "Direct sequence was refined repeatedly");
  for (double distance : {0.05, 0.12, 0.121}) {
    const Receiver nearby{{distance, 0, 0}, Isotropic{}};
    const auto result = solve(scene, tx, nearby, cfg);
    check(result.paths.size() == 1 && result.paths[0].field,
          "Nearby receiver lost its direct path");
    near(std::abs(result.impulseResponse.taps[0].coefficient), 1.0 / distance,
         "Nearby direct gain incorrect", 1e-7);
    near(result.impulseResponse.taps[0].delay, distance / C,
         "Nearby direct delay incorrect", 1e-15);
    Scene blocked;
    blocked.add(Disk{{distance / 2, 0, 0}, Vec3::UnitX(), 1}, Material::pec());
    check(solve(blocked, tx, nearby, cfg).paths.empty(),
          "Nearby receiver bypassed a blocking surface");
  }
  {
    Scene wall;
    wall.add(Rectangle{Vec3::Zero(), Vec3::UnitX(), Vec3::UnitY(), 10, 10},
             Material::pec());
    auto reflectedConfig = cfg;
    reflectedConfig.rayCount = 512;
    reflectedConfig.maxReflections = 1;
    reflectedConfig.receptionRadius = 0.5;
    const auto result = solve(wall, {{0, 0, 1}, Isotropic{}},
                              {{1, 0, 1}, Isotropic{}}, reflectedConfig);
    check(result.paths.size() == 2, "Planar reuse lost a direct or reflected path");
    std::size_t converged = 0;
    for (const auto &candidate : result.candidates)
      if (candidate.refinement.status == RefinementStatus::Converged)
        ++converged;
    check(converged == 2, "Planar reflection was refined repeatedly");
    for (const auto &ray : result.paths) {
      const auto &candidate = result.candidates.at(ray.candidate).refinement;
      near(candidate.geometry.receiver->pathDistance,
           ray.refinement.geometry.receiver->pathDistance,
           "Candidate geometry was lost when moving refinement storage");
      if (!ray.refinement.geometry.reflections.empty())
        near(ray.refinement.geometry.receiver->pathDistance, std::sqrt(5.),
             "Planar reflected path has incorrect length");
    }
  }
  tx.antenna = Isotropic{Polarization::Vertical, 2};
  const auto scaled = solve(scene, tx, rx, cfg);
  near(scaled.impulseResponse.taps[0].coefficient,
       2.0 * direct.impulseResponse.taps[0].coefficient,
       "Source directional amplitude normalized away");
  tx.antenna = ShortDipole{};
  rx.position = Vec3(0, 0, 3);
  const auto null = solve(scene, tx, rx, cfg);
  check(null.paths.size() == 1 && null.paths[0].field,
        "Source null prevented geometry discovery");
  near(null.impulseResponse.taps[0].coefficient, 0,
       "Source null generated field");

  // Colocated patch antennas must find a physical return rather than a direct ray.
  const auto pec = Material::pec();
  Scene tube;
  tube.add(Cylinder{{0, 0, -0.25}, Vec3::UnitZ(), 0.127, 0.25, false}, pec);
  tube.add(Disk{{0, 0, -0.5}, -Vec3::UnitZ(), 0.127}, pec);
  auto monostaticConfig = cfg;
  monostaticConfig.rayCount = 512;
  monostaticConfig.maxReflections = 8;
  const Transmitter monostaticTx{Vec3::Zero(), patch};
  const Receiver monostaticRx{Vec3::Zero(), patch};
  const auto returned = solve(tube, monostaticTx, monostaticRx, monostaticConfig);
  std::size_t launched = 0;
  for (auto count : returned.traceStatusCounts) launched += count;
  check(launched == monostaticConfig.rayCount,
        "Colocated solver launched a zero-length direct direction");
  check(!returned.impulseResponse.taps.empty(), "Colocated patch return missing");
  bool bottomReturn = false;
  for (const auto &ray : returned.paths) {
    check(ray.refinement.geometry.receiver->pathDistance > 0,
          "Colocated solver returned a zero-length path");
    const auto &hits = ray.refinement.geometry.reflections;
    if (hits.size() == 1 && hits[0].surfaceIndex == 1) {
      bottomReturn = true;
      near(ray.refinement.geometry.receiver->pathDistance, 1.0,
           "Bottom return should travel one metre", 1e-6);
      check(ray.field && std::isfinite(std::abs(ray.coefficient)) &&
                std::abs(ray.coefficient) > 0,
            "Colocated patch return has no finite field");
    }
  }
  check(bottomReturn, "Bottom disk reflection was not refined");
  const auto noReturn = solve(Scene{}, monostaticTx, monostaticRx, monostaticConfig);
  check(noReturn.paths.empty() && noReturn.impulseResponse.taps.empty(),
        "Empty monostatic scene generated a direct return");

  // Swapping the receive antenna reuses the incident fields without retracing.
  auto swapped = direct;
  set_receiver_antenna(swapped, ShortDipole{});
  check(swapped.paths[0].field->field == direct.paths[0].field->field,
        "Changing Rx changed the incident field");
  check(std::abs(swapped.paths[0].coefficient - direct.paths[0].coefficient) > 1e-4,
        "Rx antenna change had no effect");
  std::cout << "Antenna models, CST/HFSS imports, solver API "
               "passed\n";
  return 0;
} catch (const std::exception &e) {
  std::cerr << e.what() << '\n';
  return 1;
}
