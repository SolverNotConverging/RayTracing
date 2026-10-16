#pragma once
#include "antenna.hpp"
#include "geometry.hpp"
#include <cstdint>
#include <functional>
#include <map>
#include <string>

// Paraxial Gaussian-beam transport through material regions.
namespace rt::beam {
    using Mat2C = Eigen::Matrix2cd;
    using Frame = Eigen::Matrix<double, 3, 2>; // transverse unit axes

    // Right-handed transverse frame for a direction.
    Frame frame(const Vec3 &direction);

    enum class Termination { None, Entry, Exit };

    struct Object {
        Shape shape;
        std::optional<Material> material; // empty for terminations
        bool volume = false;
        Termination termination = Termination::None;
    };

    enum class EventKind { Interface, Terminate, Ambiguous };

    struct Event {
        double distance = 0.0;
        Vec3 position = Vec3::Zero(), normal = Vec3::UnitZ();
        Surface face;
        std::size_t object = 0;
        int before = -1, after = -1; // regions; -1 is the background
        EventKind kind = EventKind::Interface;
    };

    // Regions are the volumes; the last added volume containing a point wins.
    // Sheets (non-volume surfaces) override coincident region interfaces.
    class Scene {
    public:
        explicit Scene(Material background = {}, double tolerance = 1e-8);
        std::size_t add_volume(Shape shape, Material material); // Box, Sphere or capped Cylinder
        std::size_t add_surface(Shape shape, Material material); // conductor or PEC sheet
        // Absorb beams: surfaces on any hit, volumes on entry or exit.
        std::size_t add_termination(Shape shape, bool volume = false, Termination on = Termination::Entry);

        int region_at(const Vec3 &point) const;
        const Material &material(int region) const;
        bool starts_terminated(const Vec3 &point) const;
        std::optional<Event> next_event(const Vec3 &origin, const Vec3 &direction, double maxDistance) const;

        const Material &background() const { return background_; }
        double tolerance() const { return tolerance_; }
        const std::vector<Object> &objects() const { return objects_; }

    private:
        std::size_t add(Object object);
        Material background_;
        double tolerance_;
        std::vector<Object> objects_;
        std::vector<std::pair<std::size_t, Surface>> faces_;
    };

    enum class InteractionKind { Reflect, Transmit };

    struct Interaction {
        std::size_t object;
        InteractionKind kind;
        int before, after;
    };

    struct State {
        Vec3 position = Vec3::Zero(), direction = Vec3::UnitZ();
        Frame basis = frame(Vec3::UnitZ());
        Mat2C curvature = Mat2C::Zero(); // symmetric inverse-q matrix (1/m) in basis coordinates
        Vec3C field = Vec3C::Zero(); // on-axis E (V/m)
        int region = -1;
        double power = 0.0, distance = 0.0, opticalPath = 0.0;
        std::vector<Interaction> history;
        std::size_t branch = 0;
        std::int64_t parent = -1;

        void validate() const;
    };

    State propagate(const State &state, double distance, const Material &material, double frequencyHz);

    struct Segment {
        State state;
        double length = 0.0;
        Material material;
        double frequencyHz = 77e9;

        // Paraxial field on the half-open interval [0, length) along the axis.
        Vec3C field_at(const Vec3 &point) const;
        Eigen::Vector2d radii(double distance = 0.0) const; // 1/e field radii
    };

    struct GaussianBeam {
        Vec3 waistPosition = Vec3::Zero(), direction = Vec3::UnitZ();
        Eigen::Vector2d waistRadius = Eigen::Vector2d::Constant(0.1); // 1/e field radii
        double power = 1.0;
        Vec3C polarization = Vec3C::UnitX();
        double launchDistance = 0.0, phase = 0.0;
        std::optional<Vec3> transverseAxis;

        State launch(const Scene &scene, double frequencyHz) const;
        // Coherent decomposition into narrower beamlets (fields superpose, not powers).
        std::vector<GaussianBeam> beamlets(int samples = 17, double radiusRatio = 0.5, double extent = 4.0) const;
    };

    struct PlaneWave {
        Vec3 position = Vec3::Zero(), direction = Vec3::UnitZ();
        Vec3C polarization = Vec3C::UnitX();
        Complex amplitude = 1.0;

        State launch(const Scene &scene, double frequencyHz) const;
    };

    // Least-squares Gaussian-beamlet fit of a sampled transverse field on a plane.
    struct ApertureSource {
        Vec3 position = Vec3::Zero(), direction = Vec3::UnitZ();
        Eigen::Vector2d halfWidths = Eigen::Vector2d::Constant(0.2);
        int samples = 9;
        double beamletRadius = 0.05;
        std::function<Vec3C(const Vec3 &)> field;
        double rcond = 1e-10;
        std::optional<double> frequencyHz; // set when the field is only valid at one frequency

        std::vector<State> launch(const Scene &scene, double frequencyHz) const;
        // Sample an antenna's vacuum far field on a plane in its far zone.
        static ApertureSource from_antenna(const Antenna &antenna, const Vec3 &antennaPosition, const Vec3 &position,
                                           const Vec3 &direction, const Eigen::Vector2d &halfWidths, int samples,
                                           double beamletRadius, double frequencyHz);
    };

    using Source = std::variant<GaussianBeam, PlaneWave, ApertureSource>;

    std::vector<State> launch(const Scene &scene, const Source &source, double frequencyHz);

    struct SolverConfig {
        double frequencyHz = 77e9;
        int maxInteractions = 16;
        double maxDistance = 20.0;
        double minPowerFraction = 1e-10;
        std::size_t maxBranches = 100000; // total, including the launched states
        double maxLossRatio = 0.01; // Im(n)/Re(n) allowed for refracting dielectrics
    };

    enum class TerminalReason {
        Terminated, PowerCutoff, DistanceLimit, InteractionLimit, BranchLimit, Ambiguous, Grazing, LossyCritical
    };
    const char *to_string(TerminalReason reason);

    struct Terminal {
        std::size_t branch = 0;
        TerminalReason reason = TerminalReason::DistanceLimit;
        Vec3 position = Vec3::Zero();
        double power = 0.0;
        std::vector<Interaction> history;
    };

    struct Simulation {
        SolverConfig config;
        Scene scene;
        std::vector<Segment> segments;
        std::vector<Terminal> terminals;
        std::vector<std::string> warnings;
        double launchedPower = 0.0, absorbedPower = 0.0, interfaceFluxResidual = 0.0;

        Vec3C field_at(const Vec3 &point) const;
        std::vector<Vec3C> fields_at(const std::vector<Vec3> &points) const;
        // Constituent branch powers by fate; coherent totals need a field flux integral.
        std::map<std::string, double> power_balance() const;
    };

    // Reflected and transmitted children at an interface, or the reason the branch stops.
    std::vector<State> interact(const State &state, const Event &event, const Scene &scene, double frequencyHz,
                                std::optional<TerminalReason> &stop);

    Simulation solve_states(const Scene &scene, std::vector<State> states, const SolverConfig &config = {});
    Simulation solve(const Scene &scene, const std::vector<Source> &sources, const SolverConfig &config = {});
    Simulation solve(const Scene &scene, const Source &source, const SolverConfig &config = {});
} // namespace rt::beam
