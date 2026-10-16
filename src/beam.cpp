#include <raytracing/beam.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace rt::beam {
    namespace {
        void require(bool condition, const char *message) {
            if (!condition) throw std::invalid_argument(message);
        }

        bool positive(double x) { return std::isfinite(x) && x > 0; }

        bool transparent(const Material &m) { return m.kind == MaterialKind::Dielectric; }

        Vec3 unit(const Vec3 &v) {
            require(v.allFinite() && v.norm() > 0, "Expected a finite nonzero direction");
            return v.normalized();
        }

        Vec3C transverse_polarization(const Vec3C &p, const Vec3 &d) {
            require(p.allFinite() && p.norm() > 0, "Polarization must be a finite nonzero vector");
            require(std::abs(d.cast<Complex>().dot(p)) <= 1e-10 * p.norm(), "Polarization must be transverse");
            return p / p.norm();
        }

        Frame oriented_frame(const Vec3 &d, const std::optional<Vec3> &transverseAxis) {
            if (!transverseAxis) return frame(d);
            const Vec3 u = unit(*transverseAxis);
            require(std::abs(u.dot(d)) <= 1e-10, "Transverse axis must be perpendicular to the beam axis");
            Frame basis;
            basis << u, d.cross(u);
            return basis;
        }

        double admittance(const Material &m, double f) { return (m.index(f) / m.mu).real() / freeSpaceImpedance; }
    } // namespace

    Frame frame(const Vec3 &direction) {
        const Vec3 d = direction.normalized();
        Eigen::Index i;
        d.cwiseAbs().minCoeff(&i);
        const Vec3 u = Vec3::Unit(i).cross(d).normalized();
        Frame result;
        result << u, d.cross(u);
        return result;
    }

    // ---- Scene ---------------------------------------------------------------

    Scene::Scene(Material background, double tolerance) : background_(std::move(background)), tolerance_(tolerance) {
        background_.validate();
        require(transparent(background_), "The background must be a dielectric");
        require(positive(tolerance), "Tolerance must be positive");
    }

    std::size_t Scene::add(Object object) {
        validate_shape(object.shape);
        if (object.material) object.material->validate();
        const auto id = objects_.size();
        for (const auto &face: faces(object.shape)) faces_.emplace_back(id, face);
        objects_.push_back(std::move(object));
        return id;
    }

    std::size_t Scene::add_volume(Shape shape, Material material) {
        const auto *cylinder = std::get_if<Cylinder>(&shape);
        require(std::holds_alternative<Box>(shape) || std::holds_alternative<Sphere>(shape) ||
                (cylinder && cylinder->capped), "Volumes must be boxes, spheres or capped cylinders");
        return add({std::move(shape), std::move(material), true, Termination::None});
    }

    std::size_t Scene::add_surface(Shape shape, Material material) {
        require(!transparent(material), "Dielectrics need a volume; surfaces must be conductors or PEC");
        return add({std::move(shape), std::move(material), false, Termination::None});
    }

    std::size_t Scene::add_termination(Shape shape, bool volume, Termination on) {
        require(on != Termination::None && (volume || on == Termination::Entry),
                "Surfaces terminate on entry; volumes on entry or exit");
        return add({std::move(shape), std::nullopt, volume, on});
    }

    int Scene::region_at(const Vec3 &p) const {
        for (std::size_t i = objects_.size(); i > 0; --i) {
            const auto &o = objects_[i - 1];
            if (o.volume && o.material && contains(o.shape, p)) return static_cast<int>(i - 1);
        }
        return -1;
    }

    const Material &Scene::material(int region) const {
        return region == -1 ? background_ : *objects_.at(region).material;
    }

    bool Scene::starts_terminated(const Vec3 &p) const {
        return std::any_of(objects_.begin(), objects_.end(), [&](const auto &o) {
            return o.volume && o.termination == Termination::Entry && contains(o.shape, p);
        });
    }

    std::optional<Event> Scene::next_event(const Vec3 &origin, const Vec3 &direction, double maxDistance) const {
        struct FaceHit {
            double distance;
            std::size_t object;
            Surface face;
            Vec3 normal;
        };
        Vec3 current = origin;
        double travelled = 0;
        // Step over crossings between identical materials until something happens.
        while (maxDistance - travelled > tolerance_) {
            std::vector<FaceHit> hits;
            for (const auto &[id, face]: faces_)
                if (const auto hit = intersect(current, direction, face, tolerance_, maxDistance - travelled))
                    hits.push_back({hit->distance, id, face, hit->normal});
            if (hits.empty()) return std::nullopt;
            std::stable_sort(hits.begin(), hits.end(), [](const auto &a, const auto &b) {
                return a.distance < b.distance;
            });
            const double distance = hits[0].distance, probe = 2 * tolerance_;
            std::erase_if(hits, [&](const auto &h) { return h.distance - distance > tolerance_; });
            const Vec3 p = current + distance * direction;
            const int before = region_at(p - probe * direction), after = region_at(p + probe * direction);
            const auto event = [&](const FaceHit &h, EventKind kind) {
                return Event{travelled + distance, p, h.normal, h.face, h.object, before, after, kind};
            };

            for (const auto &h: hits) {
                const auto &o = objects_[h.object];
                if (o.termination == Termination::None) continue;
                const bool was = contains(o.shape, p - probe * direction), now = contains(o.shape, p + probe * direction);
                if (!o.volume || (o.termination == Termination::Entry && !was && now) ||
                    (o.termination == Termination::Exit && was && !now))
                    return event(h, EventKind::Terminate);
            }

            std::vector<FaceHit> active;
            for (const auto &h: hits) {
                const auto &o = objects_[h.object];
                const int id = static_cast<int>(h.object);
                if (o.termination == Termination::None && (!o.volume || id == before || id == after)) active.push_back(h);
            }
            if (!active.empty()) {
                // Sheets beat volume faces; otherwise the last added object wins.
                const FaceHit *winner = &active[0];
                bool sheet = false;
                for (const auto &h: active) {
                    const bool isSheet = !objects_[h.object].volume;
                    if ((isSheet && !sheet) || (isSheet == sheet && h.object > winner->object)) winner = &h;
                    sheet = sheet || isSheet;
                }
                for (const auto &h: active)
                    if (std::abs(winner->normal.dot(h.normal)) < 1 - 1e-10 || edge_distance(h.face, p) <= tolerance_)
                        return event(*winner, EventKind::Ambiguous);
                if (sheet || material(before) != material(after)) return event(*winner, EventKind::Interface);
            }
            current = p + probe * direction;
            travelled += distance + probe;
        }
        return std::nullopt;
    }

    // ---- States and propagation ----------------------------------------------

    void State::validate() const {
        require(position.allFinite() && std::abs(direction.norm() - 1) < 1e-10 &&
                (basis.transpose() * basis - Eigen::Matrix2d::Identity()).norm() < 1e-10 &&
                (basis.transpose() * direction).norm() < 1e-10,
                "Beam state needs a finite position, unit direction and orthonormal transverse basis");
        require(curvature.allFinite() && field.allFinite() && std::isfinite(power) && power >= 0 &&
                (curvature - curvature.transpose()).norm() <= 1e-10 * std::max(1.0, curvature.norm()),
                "Beam state needs a finite symmetric curvature, finite field and nonnegative power");
    }

    State propagate(const State &state, double distance, const Material &material, double frequencyHz) {
        const Complex n = material.index(frequencyHz);
        const double k0 = wavenumber(frequencyHz);
        const auto &h = state.curvature;
        // det(I + distance h) through the eigenvalues of h, with a continuous square root.
        const Complex halfTrace = h.trace() / 2.0;
        const Complex delta = std::sqrt((h(0, 0) - h(1, 1)) * (h(0, 0) - h(1, 1)) / 4.0 + h(0, 1) * h(1, 0));
        const Complex a = 1.0 + distance * (halfTrace + delta), b = 1.0 + distance * (halfTrace - delta);
        require(std::min(std::abs(a), std::abs(b)) >= 1e-14,
                "Plane wave reached a geometrical focus; use a finite GaussianBeam");
        State result = state;
        result.position += distance * state.direction;
        result.curvature = h * (Mat2C::Identity() + distance * h).inverse();
        result.field = state.field / (std::sqrt(a) * std::sqrt(b)) * std::exp(Complex(0, 1) * k0 * n * distance);
        result.power *= std::exp(-2 * k0 * n.imag() * distance);
        result.distance += distance;
        result.opticalPath += n.real() * distance;
        return result;
    }

    Vec3C Segment::field_at(const Vec3 &point) const {
        const Vec3 offset = point - state.position;
        const double s = offset.dot(state.direction);
        if (s < 0 || s >= length) return Vec3C::Zero();
        const auto at = propagate(state, s, material, frequencyHz);
        const Eigen::Vector2cd x = (state.basis.transpose() * offset).cast<Complex>();
        const double k = wavenumber(frequencyHz) * material.index(frequencyHz).real();
        return at.field * std::exp(Complex(0, 0.5) * k * (x.transpose() * at.curvature * x).value());
    }

    Eigen::Vector2d Segment::radii(double distance) const {
        const auto at = propagate(state, distance, material, frequencyHz);
        const Eigen::Vector2d values = Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>(at.curvature.imag()).eigenvalues();
        if (values.minCoeff() <= 0) return Eigen::Vector2d::Constant(infinity);
        const double k = wavenumber(frequencyHz) * material.index(frequencyHz).real();
        return (2 / (k * values.array())).sqrt();
    }

    // ---- Sources -------------------------------------------------------------

    State GaussianBeam::launch(const Scene &scene, double frequencyHz) const {
        const Vec3 d = unit(direction);
        const Frame basis = oriented_frame(d, transverseAxis);
        require(waistPosition.allFinite() && (waistRadius.array() > 0).all() && positive(power) &&
                std::isfinite(phase) && std::isfinite(launchDistance),
                "Gaussian beam needs a finite waist, positive radii and power, and finite phase");
        const Vec3 origin = waistPosition + launchDistance * d;
        const int region = scene.region_at(origin);
        const auto &medium = scene.material(region);
        const Complex n = medium.index(frequencyHz);
        require(n.imag() == 0, "Launch Gaussian beams in a lossless dielectric region");
        const double k = wavenumber(frequencyHz) * n.real();
        require((2 / (k * waistRadius.array())).maxCoeff() <= 0.2,
                "Gaussian divergence exceeds 0.2 rad, outside the paraxial model");
        const Eigen::Vector2d rayleigh = k * waistRadius.array().square() / 2;
        Mat2C h = Mat2C::Zero();
        Complex gouy = 1;
        for (int i = 0; i < 2; ++i) {
            const Complex q(launchDistance, -rayleigh[i]);
            h(i, i) = 1.0 / q;
            gouy *= std::sqrt(Complex(0, -rayleigh[i]) / q);
        }
        const double amplitude = std::sqrt(4 * power / (pi * waistRadius.prod() * admittance(medium, frequencyHz)));
        State result;
        result.position = origin;
        result.direction = d;
        result.basis = basis;
        result.curvature = h;
        result.field = amplitude * gouy * std::exp(Complex(0, k * launchDistance + phase)) *
                       transverse_polarization(polarization, d);
        result.region = region;
        result.power = power;
        return result;
    }

    std::vector<GaussianBeam> GaussianBeam::beamlets(int samples, double radiusRatio, double extent) const {
        require(samples >= 3 && radiusRatio > 0 && radiusRatio < 1 && positive(extent),
                "Beamlets need samples >= 3, 0 < radius_ratio < 1 and a positive extent");
        const Vec3 d = unit(direction);
        const Frame basis = oriented_frame(d, transverseAxis);
        const Eigen::Vector2d small = waistRadius * radiusRatio;
        const Eigen::Vector2d spread = (waistRadius.array().square() - small.array().square()).sqrt();
        const Eigen::Vector2d step = 2 * extent * spread / (samples - 1);
        const auto trapezoid = [&](int i) { return i == 0 || i == samples - 1 ? 0.5 : 1.0; };
        std::vector<GaussianBeam> result;
        for (int i = 0; i < samples; ++i)
            for (int j = 0; j < samples; ++j) {
                const Eigen::Vector2d xy(-extent * spread[0] + i * step[0], -extent * spread[1] + j * step[1]);
                const double weight = step.prod() * trapezoid(i) * trapezoid(j);
                const double amplitude = weight * waistRadius.prod() / (pi * small.prod() * spread.prod()) *
                                         std::exp(-(xy.array() / spread.array()).square().sum());
                GaussianBeam beam = *this;
                beam.waistPosition += basis * xy;
                beam.direction = d;
                beam.waistRadius = small;
                beam.power = power * amplitude * amplitude * small.prod() / waistRadius.prod();
                beam.transverseAxis = basis.col(0);
                result.push_back(std::move(beam));
            }
        return result;
    }

    State PlaneWave::launch(const Scene &scene, double frequencyHz) const {
        const Vec3 d = unit(direction);
        const int region = scene.region_at(position);
        const auto &medium = scene.material(region);
        State result;
        result.position = position;
        result.direction = d;
        result.basis = frame(d);
        result.field = amplitude * transverse_polarization(polarization, d);
        result.region = region;
        result.power = std::norm(amplitude) * admittance(medium, frequencyHz) / 2;
        return result;
    }

    std::vector<State> ApertureSource::launch(const Scene &scene, double f) const {
        require(!frequencyHz || f == *frequencyHz, "This aperture field was sampled at a different frequency");
        require(samples >= 3 && (halfWidths.array() > 0).all() && positive(beamletRadius) && rcond > 0 && rcond < 1 &&
                field, "Aperture needs samples >= 3, positive sizes, 0 < rcond < 1 and a field callback");
        const Vec3 d = unit(direction);
        const Frame basis = frame(d);
        const int count = samples * samples;
        Eigen::MatrixXd xy(count, 2), real(count, 3), imag(count, 3);
        std::vector<Vec3> positions;
        for (int i = 0; i < samples; ++i)
            for (int j = 0; j < samples; ++j) {
                const int row = i * samples + j;
                xy.row(row) << -halfWidths[0] + 2 * halfWidths[0] * i / (samples - 1),
                        -halfWidths[1] + 2 * halfWidths[1] * j / (samples - 1);
                positions.push_back(position + basis * xy.row(row).transpose());
                const Vec3C e = field(positions.back());
                real.row(row) = e.real();
                imag.row(row) = e.imag();
            }
        // Least-squares beamlet weights: sum_j K(x_i - x_j) w_j = E(x_i).
        Eigen::MatrixXd kernel(count, count);
        for (int i = 0; i < count; ++i)
            for (int j = 0; j < count; ++j)
                kernel(i, j) = std::exp(-(xy.row(i) - xy.row(j)).squaredNorm() / (beamletRadius * beamletRadius));
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(kernel, Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd.setThreshold(rcond);
        const Eigen::MatrixXd weightsReal = svd.solve(real), weightsImag = svd.solve(imag);

        const auto &medium = scene.material(scene.region_at(position));
        require(medium.index(f).imag() == 0, "Aperture sources need a lossless launch region");
        std::vector<State> states;
        for (int i = 0; i < count; ++i) {
            const Vec3C e = weightsReal.row(i).transpose().cast<Complex>() +
                            Complex(0, 1) * weightsImag.row(i).transpose().cast<Complex>();
            const double amplitude = e.norm();
            if (amplitude <= 1e-30) continue;
            GaussianBeam beam;
            beam.waistPosition = positions[i];
            beam.direction = d;
            beam.waistRadius.setConstant(beamletRadius);
            beam.power = pi * beamletRadius * beamletRadius * admittance(medium, f) / 4 * amplitude * amplitude;
            beam.polarization = e / amplitude;
            beam.transverseAxis = basis.col(0);
            states.push_back(beam.launch(scene, f));
        }
        require(!states.empty(), "Aperture field is zero");
        return states;
    }

    ApertureSource ApertureSource::from_antenna(const Antenna &antenna, const Vec3 &antennaPosition,
                                                const Vec3 &position, const Vec3 &direction,
                                                const Eigen::Vector2d &halfWidths, int samples, double beamletRadius,
                                                double frequencyHz) {
        const Vec3 d = unit(direction);
        ApertureSource source{position, d, halfWidths, samples, beamletRadius, {}, 1e-10, frequencyHz};
        source.field = [antenna, antennaPosition, d, frequencyHz](const Vec3 &p) -> Vec3C {
            const Vec3 offset = p - antennaPosition;
            const double distance = offset.norm();
            const Vec3C e = antenna.farfield(offset / distance, frequencyHz) / distance;
            const Vec3C dc = d.cast<Complex>();
            return (e - dc * dc.dot(e)) * std::exp(Complex(0, wavenumber(frequencyHz) * distance));
        };
        return source;
    }

    std::vector<State> launch(const Scene &scene, const Source &source, double frequencyHz) {
        return std::visit([&](const auto &s) -> std::vector<State> {
            if constexpr (std::is_same_v<std::decay_t<decltype(s)>, ApertureSource>) return s.launch(scene, frequencyHz);
            else return {s.launch(scene, frequencyHz)};
        }, source);
    }

    // ---- Interfaces ----------------------------------------------------------

    std::vector<State> interact(const State &state, const Event &event, const Scene &scene, double f,
                                std::optional<TerminalReason> &stop) {
        stop.reset();
        const Vec3 &d = state.direction;
        const Frame &u = state.basis;
        const Vec3 normal = event.normal.dot(d) > 0 ? Vec3(-event.normal) : event.normal;
        const double c = -normal.dot(d);
        if (c < 1e-10) {
            stop = TerminalReason::Grazing;
            return {};
        }
        const auto &incoming = scene.material(state.region);
        const auto &object = scene.objects()[event.object];
        const auto &outgoing = object.volume ? scene.material(event.after) : *object.material;
        const auto coefficients = fresnel(incoming, outgoing, c, f);
        // s: perpendicular to the plane of incidence; p = s x d.
        Vec3 s = d.cross(normal);
        s = s.norm() > 1e-10 ? s.normalized() : Vec3(u.col(0));
        const Vec3 p = s.cross(d);
        const Complex es = s.cast<Complex>().dot(state.field), ep = p.cast<Complex>().dot(state.field);
        // Transport the transverse ray bundle across the curved interface.
        const Eigen::Matrix3d projection = Eigen::Matrix3d::Identity() - d * normal.transpose() / normal.dot(d);
        const Frame dx = projection * u, dn = curvature(event.face, event.position, normal) * dx;

        std::vector<State> children;
        for (const auto kind: {InteractionKind::Reflect, InteractionKind::Transmit}) {
            Vec3 direction;
            Frame derivativeX, derivativeD;
            Complex cs, cp;
            int region;
            double powerRatio = 1;
            if (kind == InteractionKind::Reflect) {
                direction = d + 2 * c * normal;
                cs = coefficients.rs;
                cp = coefficients.rp;
                derivativeX = -2 * (normal * (d.transpose() * dn) + d.dot(normal) * dn);
                derivativeD = (Eigen::Matrix3d::Identity() - 2 * normal * normal.transpose()) * u;
                region = state.region;
            } else {
                if (!transparent(outgoing)) continue;
                const double n1 = incoming.index(f).real(), n2 = outgoing.index(f).real();
                const double eta = n1 / n2, ct2 = 1 - eta * eta * (1 - c * c);
                if (ct2 <= 1e-14) { // total internal reflection
                    if (outgoing.index(f).imag() != 0 || incoming.index(f).imag() != 0) {
                        stop = TerminalReason::LossyCritical;
                        return {};
                    }
                    continue;
                }
                const double ct = std::sqrt(ct2), derivative = eta - eta * eta * c / ct;
                direction = eta * d + (eta * c - ct) * normal;
                cs = coefficients.ts;
                cp = coefficients.tp;
                derivativeX = (eta * c - ct) * dn - normal * (derivative * (d.transpose() * dn));
                derivativeD = eta * u - normal * (derivative * (normal.transpose() * u));
                region = event.after;
                powerRatio = admittance(outgoing, f) * std::abs(direction.normalized().dot(normal)) /
                             (admittance(incoming, f) * c);
            }
            direction.normalize();
            const Frame v = frame(direction);
            const Eigen::Matrix2d a = v.transpose() * dx;
            Mat2C h = ((v.transpose() * derivativeX).cast<Complex>() +
                       (v.transpose() * derivativeD).cast<Complex>() * state.curvature) *
                      a.inverse().cast<Complex>();
            h = ((h + h.transpose()) / 2.0).eval(); // symmetric up to rounding
            const Vec3C field = cs * es * s.cast<Complex>() + cp * ep * s.cross(direction).cast<Complex>();
            State child = state;
            child.position = event.position;
            child.direction = direction;
            child.basis = v;
            child.curvature = h;
            child.field = field;
            child.region = region;
            child.power = state.power * powerRatio * field.squaredNorm() / std::max(state.field.squaredNorm(), 1e-300);
            child.history.push_back({event.object, kind, state.region, region});
            children.push_back(std::move(child));
        }
        return children;
    }

    // ---- Solver --------------------------------------------------------------

    const char *to_string(TerminalReason reason) {
        switch (reason) {
            case TerminalReason::Terminated: return "terminated";
            case TerminalReason::PowerCutoff: return "power_cutoff";
            case TerminalReason::DistanceLimit: return "distance_limit";
            case TerminalReason::InteractionLimit: return "interaction_limit";
            case TerminalReason::BranchLimit: return "branch_limit";
            case TerminalReason::Ambiguous: return "ambiguous";
            case TerminalReason::Grazing: return "grazing";
            case TerminalReason::LossyCritical: return "lossy_critical";
        }
        return "unknown";
    }

    Vec3C Simulation::field_at(const Vec3 &p) const {
        const auto &medium = scene.material(scene.region_at(p));
        if (!transparent(medium) || scene.starts_terminated(p)) return Vec3C::Zero();
        Vec3C result = Vec3C::Zero();
        for (const auto &segment: segments)
            if (segment.material == medium) result += segment.field_at(p);
        return result;
    }

    std::vector<Vec3C> Simulation::fields_at(const std::vector<Vec3> &points) const {
        std::vector<Vec3C> values;
        values.reserve(points.size());
        for (const auto &p: points) values.push_back(field_at(p));
        return values;
    }

    ImpulseResponse Simulation::impulse_response(const Vec3 &point, const Vec3C &analyzer,
                                                 double delayTolerance) const {
        std::vector<ImpulseTap> taps;
        const auto &medium = scene.material(scene.region_at(point));
        if (!transparent(medium) || scene.starts_terminated(point)) return {config.frequencyHz, {}};
        const double n = medium.index(config.frequencyHz).real();
        for (std::size_t i = 0; i < segments.size(); ++i) {
            const auto &segment = segments[i];
            if (segment.material != medium) continue;
            const Vec3C field = segment.field_at(point); // zero outside the segment
            if (field.isZero(0)) continue;
            const double distance = (point - segment.state.position).dot(segment.state.direction);
            const double delay = (segment.state.opticalPath + n * distance) / speedOfLight;
            taps.push_back({delay, analyzer.dot(field), field, {i}});
        }
        return group_taps(std::move(taps), config.frequencyHz, delayTolerance);
    }

    std::map<std::string, double> Simulation::power_balance() const {
        std::map<std::string, double> result{
            {"launched", launchedPower}, {"absorbed", absorbedPower}, {"interface_flux_residual", interfaceFluxResidual}};
        for (const auto &t: terminals) result[to_string(t.reason)] += t.power;
        return result;
    }

    Simulation solve_states(const Scene &scene, std::vector<State> states, const SolverConfig &config) {
        const double f = config.frequencyHz;
        require(positive(f) && positive(config.maxDistance) && config.maxInteractions >= 0 && config.maxBranches >= 1 &&
                config.minPowerFraction >= 0 && config.minPowerFraction < 1,
                "Invalid beam solver frequency, distance, budgets or power threshold");
        require(!states.empty() && states.size() <= config.maxBranches,
                "Need at least one source state, and no more than max_branches");
        const auto weakLoss = [&](const Material &m) {
            if (!transparent(m)) return;
            const Complex n = m.index(f);
            require(n.imag() <= config.maxLossRatio * n.real(),
                    "Dielectric loss exceeds max_loss_ratio; strongly absorbing refraction is unsupported");
        };
        weakLoss(scene.background());
        for (const auto &o: scene.objects())
            if (o.material) weakLoss(*o.material);

        Simulation result{config, scene};
        for (std::size_t i = 0; i < states.size(); ++i) {
            states[i].validate();
            require(transparent(scene.material(states[i].region)), "Cannot launch inside an opaque material");
            states[i].branch = i;
            result.launchedPower += states[i].power;
        }
        const auto warn = [&](const std::string &message) {
            if (std::find(result.warnings.begin(), result.warnings.end(), message) == result.warnings.end())
                result.warnings.push_back(message);
        };
        if (states.size() > 1)
            warn("Coherent sources: power_balance sums branch powers; integrate the field flux for physical totals");
        const auto finish = [&](const State &s, TerminalReason reason) {
            result.terminals.push_back({s.branch, reason, s.position, s.power, s.history});
        };

        std::size_t allocated = states.size(), nextId = allocated;
        auto &queue = states;
        while (!queue.empty()) {
            State state = std::move(queue.back());
            queue.pop_back();
            if (scene.starts_terminated(state.position)) {
                finish(state, TerminalReason::Terminated);
                continue;
            }
            if (state.power == 0 || state.power < config.minPowerFraction * result.launchedPower) {
                finish(state, TerminalReason::PowerCutoff);
                continue;
            }
            const double remaining = config.maxDistance - state.distance;
            if (remaining <= scene.tolerance()) {
                finish(state, TerminalReason::DistanceLimit);
                continue;
            }
            const auto event = scene.next_event(state.position, state.direction, remaining);
            const double length = event ? event->distance : remaining;
            const auto &material = scene.material(state.region);
            const Segment segment{state, length, material, f};
            result.segments.push_back(segment);
            State end = propagate(state, length, material, f);
            result.absorbedPower += state.power - end.power;
            if (!event) {
                finish(end, TerminalReason::DistanceLimit);
                continue;
            }
            if (event->kind != EventKind::Interface) {
                finish(end, event->kind == EventKind::Ambiguous ? TerminalReason::Ambiguous : TerminalReason::Terminated);
                continue;
            }
            if (state.history.size() >= static_cast<std::size_t>(config.maxInteractions)) {
                finish(end, TerminalReason::InteractionLimit);
                continue;
            }
            end.region = event->before;
            if (const auto widths = segment.radii(length); widths.allFinite()) {
                const double footprint = 3 * widths.maxCoeff() / std::max(std::abs(state.direction.dot(event->normal)), 1e-10);
                if (edge_distance(event->face, event->position) < footprint)
                    warn("Branch " + std::to_string(state.branch) + ": beam footprint approaches object " +
                         std::to_string(event->object) + " edge; clipping is not modelled");
            }
            std::optional<TerminalReason> stop;
            auto children = interact(end, *event, scene, f, stop);
            if (stop) {
                finish(end, *stop);
                continue;
            }
            double outgoingPower = 0;
            for (const auto &child: children) outgoingPower += child.power;
            const double residual = end.power - outgoingPower;
            const auto &target = scene.objects()[event->object];
            const auto &boundary = target.volume ? scene.material(event->after) : *target.material;
            if (!transparent(boundary)) result.absorbedPower += residual;
            else {
                result.interfaceFluxResidual += residual;
                if (std::abs(residual) > 1e-7 * std::max(end.power, 1e-300))
                    warn("Lossy-interface branch budgets omit interference flux; see interface_flux_residual");
            }
            for (auto &child: children) {
                child.parent = static_cast<std::int64_t>(state.branch);
                child.branch = nextId++;
                if (allocated >= config.maxBranches) finish(child, TerminalReason::BranchLimit);
                else {
                    ++allocated;
                    queue.push_back(std::move(child));
                }
            }
        }
        return result;
    }

    Simulation solve(const Scene &scene, const std::vector<Source> &sources, const SolverConfig &config) {
        std::vector<State> states;
        for (const auto &source: sources) {
            auto launched = launch(scene, source, config.frequencyHz);
            states.insert(states.end(), launched.begin(), launched.end());
        }
        return solve_states(scene, std::move(states), config);
    }

    Simulation solve(const Scene &scene, const Source &source, const SolverConfig &config) {
        return solve(scene, std::vector<Source>{source}, config);
    }
} // namespace rt::beam
