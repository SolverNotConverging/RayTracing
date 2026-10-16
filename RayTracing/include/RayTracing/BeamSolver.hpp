#pragma once
#include "Antenna.hpp"
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <tuple>

// Paraxial finite-beam transport, SI units, exp(-i omega t).
// Independent of Python, persistence formats and graphics backends.
namespace rt::beam {
using Mat2C = Eigen::Matrix2cd;
using Frame = Eigen::Matrix<double, 3, 2>;
// Passive isotropic material. Positive Im(epsilon) and conductivity represent loss.
struct Material {
    std::string kind = "dielectric";
    Complex epsilon = 1;
    double mu = 1, conductivity = 0;
    std::string name = "vacuum";
    void validate() const;
    Complex permittivity(double frequency) const;
    Complex index(double frequency) const;
    bool operator==(const Material &) const = default;
};
std::array<Complex, 4> fresnel(const Material &, const Material &, double cosine, double frequency);
struct Box {
    Vec3 center = Vec3::Zero(), halfLengths = Vec3::Ones();
    Eigen::Matrix3d axes = Eigen::Matrix3d::Identity();
    void validate() const;
};
using Shape = std::variant<Box, Rectangle, Disk, Sphere, Cylinder, Triangle>;
std::vector<Surface> faces(const Shape &);
bool contains(const Shape &, const Vec3 &);
double edge_distance(const Surface &, const Vec3 &);
Eigen::Matrix3d curvature(const Surface &, const Vec3 &, const Vec3 &normal);
Vec3 unit(const Vec3 &);
Frame frame(const Vec3 &);
Vec3C polarization(const Vec3C &, const Vec3 &direction);
struct Object {
    Shape shape;
    std::optional<Material> material;
    bool volume = false;
    std::string terminate;
};
struct Event {
    double distance = 0;
    Vec3 position = Vec3::Zero(), normal = Vec3::UnitZ();
    Surface shape;
    std::size_t objectId = 0;
    int before = -1, after = -1;
    std::string kind;
};
// Last-added containing volume wins. Opaque sheets override coincident region interfaces.
// Geometry/materials are validated and copied at insertion; accessors are const.
class Scene {
public:
    explicit Scene(Material background = {}, double tolerance = 1e-8);
    std::size_t add(Shape, std::optional<Material>, bool volume, std::string terminate = "");
    int region_at(const Vec3 &) const;
    const Material &material(int region) const;
    bool starts_terminated(const Vec3 &) const;
    std::optional<Event> next_event(const Vec3 &, const Vec3 &direction, double maxDistance) const;
    void validate() const;
    const Material &background() const { return background_; }
    double tolerance() const { return tolerance_; }
    const std::vector<Object> &objects() const { return objects_; }
private:
    Material background_;
    double tolerance_;
    std::vector<Object> objects_;
    std::vector<std::pair<std::size_t, Surface>> faces_;
};
// (object index, reflect/transmit, incident region, outgoing region); -1 is background.
using Interaction = std::tuple<std::size_t, std::string, int, int>;
struct State {
    Vec3 position = Vec3::Zero(), direction = Vec3::UnitZ();
    Frame basis = frame(Vec3::UnitZ());
    // Symmetric inverse-q matrix (1/m), expressed in the transverse basis columns.
    Mat2C curvature = Mat2C::Zero();
    Vec3C field = Vec3::Zero().cast<Complex>();
    int region = -1;
    // Constituent power (W), travelled distance (m), integral of Re(n) ds (m).
    double power = 0, distance = 0, opticalPath = 0;
    std::vector<Interaction> history;
    std::size_t branch = 0;
    std::int64_t parent = -1;
    void validate() const;
};
State propagate(const State &, double distance, const Material &, double frequency);
struct Segment {
    State state;
    double length = 0;
    Material material;
    double frequency = 77e9;
    // Complex Cartesian E (V/m); segment intervals are half open to avoid double counting.
    Vec3C field_at(const Vec3 &) const;
    Eigen::Vector2d radii(double distance = 0) const;
};
struct GaussianBeam {
    Vec3 waistPosition = Vec3::Zero(), direction = Vec3::UnitZ();
    // 1/e field radii along the transverse frame axes, in metres.
    Eigen::Vector2d waistRadius = Eigen::Vector2d::Constant(.1);
    double power = 1;
    Vec3C polarization = Vec3::UnitX().cast<Complex>();
    double launchDistance = 0, phase = 0;
    std::optional<Vec3> transverseAxis;
    State launch(const Scene &, double frequency) const;
    std::vector<GaussianBeam> beamlets(int samples = 17, double radiusRatio = .5, double extent = 4) const;
};
struct PlaneWave {
    Vec3 position = Vec3::Zero(), direction = Vec3::UnitZ();
    Vec3C polarization = Vec3::UnitX().cast<Complex>();
    Complex amplitude = 1;
    State launch(const Scene &, double frequency) const;
};
struct ApertureSource {
    Vec3 position = Vec3::Zero(), direction = Vec3::UnitZ();
    Eigen::Vector2d halfWidths = Eigen::Vector2d::Constant(.2);
    int samples = 9;
    double beamletRadius = .05;
    std::function<Vec3C(const Vec3 &)> field;
    double rcond = 1e-10;
    std::optional<double> frequencyHz;
    bool vacuumOnly = false;
    std::vector<State> launch(const Scene &, double frequency) const;
    static ApertureSource from_antenna(const rt::Antenna &, const Vec3 &antennaPosition,
        const Vec3 &position, const Vec3 &direction, const Eigen::Vector2d &halfWidths,
        int samples, double beamletRadius, double frequency);
};
using Source = std::variant<GaussianBeam, PlaneWave, ApertureSource>;
struct SolverConfig {
    double frequencyHz = 77e9;
    int maxInteractions = 16;
    double maxDistance = 20, minPowerFraction = 1e-10;
    // Total allocated branches, including initial sources, rather than queue length.
    std::size_t maxBranches = 100000;
    double maxLossRatio = .01;
    void validate() const;
};
struct Terminal {
    std::size_t branch = 0;
    std::string reason;
    Vec3 position = Vec3::Zero();
    double power = 0;
    std::vector<Interaction> history;
};
struct Simulation {
    SolverConfig config;
    Scene scene;
    std::vector<Segment> segments;
    std::vector<Terminal> terminals;
    std::vector<std::string> warnings;
    double launchedPower = 0, absorbedPower = 0, interfaceFluxResidual = 0;
    // Complex Cartesian E (V/m); segment intervals are half open to avoid double counting.
    Vec3C field_at(const Vec3 &) const;
    std::vector<Vec3C> fields_at(const std::vector<Vec3> &) const;
    // Constituent branch powers; coherent total power requires a field-flux integral.
    std::map<std::string, double> power_balance() const;
};
std::pair<std::vector<State>, std::string> interact(const State &, const Event &, const Scene &, double frequency);
// solve_states supports custom sources producing validated launch states.
Simulation solve_states(const Scene &, const std::vector<State> &, const SolverConfig & = {});
Simulation solve(const Scene &, const std::vector<Source> &, const SolverConfig & = {});
Simulation solve(const Scene &, const Source &, const SolverConfig & = {});
} // namespace rt::beam
