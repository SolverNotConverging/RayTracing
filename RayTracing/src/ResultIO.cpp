#include "ResultIO.hpp"
#include "Constants.hpp"
#include <algorithm>
#include <fstream>
#include <hdf5.h>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <sstream>
#include <stdexcept>

namespace nlohmann {
    template<>
    struct adl_serializer<Complex> {
        static void to_json(json &j, const Complex &c) {
            j = json::array({c.real(), c.imag()});
        }

        static void from_json(const json &j, Complex &c) {
            c = {j.at(0).get<double>(), j.at(1).get<double>()};
        }
    };

    template<class T>
    struct adl_serializer<std::optional<T> > {
        static void to_json(json &j, const std::optional<T> &v) {
            if (v)
                j = *v;
            else
                j = nullptr;
        }

        static void from_json(const json &j, std::optional<T> &v) {
            if (j.is_null())
                v.reset();
            else
                v = j.get<T>();
        }
    };

    template<class T, int R, int C, int O, int MR, int MC>
    struct adl_serializer<Eigen::Matrix<T, R, C, O, MR, MC> > {
        using M = Eigen::Matrix<T, R, C, O, MR, MC>;

        static void to_json(json &j, const M &m) {
            j = json::array();
            for (int r = 0; r < R; ++r)
                for (int c = 0; c < C; ++c)
                    j.push_back(m(r, c));
        }

        static void from_json(const json &j, M &m) {
            if (!j.is_array() || j.size() != R * C)
                throw std::runtime_error("Invalid saved matrix dimensions");
            std::size_t i = 0;
            for (int r = 0; r < R; ++r)
                for (int c = 0; c < C; ++c)
                    m(r, c) = j.at(i++).get<T>();
        }
    };
} // namespace nlohmann
using nlohmann::json;

namespace rt {
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Material, name, perfectConductor, relativePermittivity,
                                       lossTangent, conductivity, relativePermeability)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ReflectionCoefficients, incidenceAngle, s, p)
} // namespace rt

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Rectangle, center_, u_, v_, halfWidth_,
                                   halfHeight_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Sphere, center_, radius_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Disk, center_, normal_, radius_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Cylinder, center_, axis_, radius_,
                                   halfLength_, capped_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Triangle, a_, b_, c_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SequenceReflection, surfaceIndex_, position_,
                                   normal_, segmentDistance_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ReceiverApproach, closestPoint_, residual_,
                                   missDistance_, finalSegmentDistance_,
                                   pathDistance_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SequenceEvaluation, status_, reflections_,
                                   finalOrigin_, finalDirection_, receiver_,
                                   blockingSurfaceIndex_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RefinementResult, status_,
                                   transmitterPosition_, receiverPosition_,
                                   launchDirection_, geometry_, iterations_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SpreadingOptions, angularStep_,
                                   maxStepHalvings_, receiverTolerance_,
                                   relativeDerivativeTolerance_,
                                   absoluteDerivativeTolerance_,
                                   minimumSingularValue_,
                                   minimumSingularValueRatio_,
                                   referenceDistance_, tMin_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SpreadingResult, status_, jacobian_,
                                   singularValues_, areaPerSolidAngle_,
                                   angularStep_, derivativeDifference_,
                                   fieldFactor_, launchU_, launchV_, receiverU_,
                                   receiverV_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(DeduplicationOptions, positionTolerance_,
                                   angleTolerance_, lengthTolerance_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ReconstructedField, receiverField_,
                                   normalizedReceiverField_,
                                   transportedReferenceField_,
                                   arrivalDirection_, reflections_, pathDistance_,
                                   opticalPath_, delaySeconds_, frequencyHz_,
                                   fieldFactor_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ImpulseTap, delaySeconds_, coefficient_,
                                   normalizedField_, pathIndices_)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ImpulseResponse, frequencyHz_,
                                   receiverPolarization_, taps_,
                                   delayToleranceSeconds_, receiverModel_)

void to_json(json &j, const RefinementOptions &o) {
    j = {
        {"maxIterations_", o.maxIterations_},
        {"receiverTolerance_", o.receiverTolerance_},
        {"finiteDifferenceStep_", o.finiteDifferenceStep_},
        {"maxDirectionStep_", o.maxDirectionStep_},
        {"maxBacktracks_", o.maxBacktracks_},
        {"tMin_", o.tMin_},
        {"cornerSeparationTolerance_", o.cornerSeparationTolerance_},
        {
            "maxPathDistance_", std::isfinite(o.maxPathDistance_)
                                    ? json(o.maxPathDistance_)
                                    : json(nullptr)
        }
    };
}

void from_json(const json &j, RefinementOptions &o) {
    j.at("maxIterations_").get_to(o.maxIterations_);
    j.at("receiverTolerance_").get_to(o.receiverTolerance_);
    j.at("finiteDifferenceStep_").get_to(o.finiteDifferenceStep_);
    j.at("maxDirectionStep_").get_to(o.maxDirectionStep_);
    j.at("maxBacktracks_").get_to(o.maxBacktracks_);
    j.at("tMin_").get_to(o.tMin_);
    j.at("cornerSeparationTolerance_").get_to(o.cornerSeparationTolerance_);
    o.maxPathDistance_ = j.at("maxPathDistance_").is_null()
                             ? std::numeric_limits<double>::infinity()
                             : j.at("maxPathDistance_").get<double>();
}

namespace rt {
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Medium, relativePermittivity,
                                       relativePermeability)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Isotropic, polarization, amplitude)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ShortDipole, effectiveLengthMetres,
                                       currentAmperes)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ThinWireDipole, lengthMetres,
                                       feedCurrentAmperes)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RectangularAperture, widthMetres,
                                       heightMetres, fieldX, fieldY, pecBacked)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SolverConfig, frequencyHz, rayCount,
                                       maxReflections, maxDistance, receptionRadius,
                                       commonSourceReference, delayToleranceSeconds,
                                       medium, refinement, spreading, deduplication)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Candidate, refinement, coarseVertices)
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SolvedRay, candidateIndex, refinement,
                                       spreading, field, sourceField,
                                       receivedCoefficient)

    namespace {
        // Version 2: per-surface materials, reflection coefficients, exact c0.
        constexpr int schemaVersion = 2;

        struct H5 {
            hid_t id = -1;

            herr_t (*close)(hid_t) = nullptr;

            H5(hid_t value, herr_t (*closer)(hid_t)) : id(value), close(closer) {
                if (id < 0)
                    throw std::runtime_error("HDF5 operation failed");
            }

            ~H5() {
                if (id >= 0)
                    close(id);
            }

            H5(const H5 &) = delete;

            H5 &operator=(const H5 &) = delete;
        };

        void ensure(herr_t status) {
            if (status < 0)
                throw std::runtime_error("HDF5 read/write failed");
        }

        void group(hid_t file, const std::string &name) {
            H5 g(H5Gcreate2(file, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                 H5Gclose);
        }

        void write_text(hid_t file, const std::string &name, const json &value) {
            const auto text = value.dump();
            H5 type(H5Tcopy(H5T_C_S1), H5Tclose);
            ensure(H5Tset_size(type.id, std::max<std::size_t>(1, text.size())));
            ensure(H5Tset_cset(type.id, H5T_CSET_UTF8));
            H5 space(H5Screate(H5S_SCALAR), H5Sclose);
            H5 data(H5Dcreate2(file, name.c_str(), type.id, space.id, H5P_DEFAULT,
                               H5P_DEFAULT, H5P_DEFAULT),
                    H5Dclose);
            ensure(
                H5Dwrite(data.id, type.id, H5S_ALL, H5S_ALL, H5P_DEFAULT, text.data()));
        }

        json read_text(hid_t file, const std::string &name) {
            H5 data(H5Dopen2(file, name.c_str(), H5P_DEFAULT), H5Dclose);
            H5 type(H5Dget_type(data.id), H5Tclose);
            H5 space(H5Dget_space(data.id), H5Sclose);
            if (H5Tget_class(type.id) != H5T_STRING || H5Tis_variable_str(type.id) > 0 ||
                H5Sget_simple_extent_ndims(space.id) != 0)
                throw std::runtime_error("Unexpected HDF5 metadata type");
            const std::size_t n = H5Tget_size(type.id);
            if (n > 64 * 1024 * 1024)
                throw std::runtime_error("HDF5 metadata too large");
            std::string text(n, '\0');
            ensure(H5Dread(data.id, type.id, H5S_ALL, H5S_ALL, H5P_DEFAULT, text.data()));
            return json::parse(text);
        }

        void write_numbers(hid_t file, const std::string &name,
                           const std::vector<double> &values,
                           const std::vector<hsize_t> &dims) {
            H5 space(
                H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr),
                H5Sclose);
            H5 data(H5Dcreate2(file, name.c_str(), H5T_IEEE_F64LE, space.id, H5P_DEFAULT,
                               H5P_DEFAULT, H5P_DEFAULT),
                    H5Dclose);
            if (!values.empty())
                ensure(H5Dwrite(data.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                                values.data()));
        }

        std::vector<double> read_numbers(hid_t file, const std::string &name,
                                         const std::vector<hsize_t> &expected) {
            H5 data(H5Dopen2(file, name.c_str(), H5P_DEFAULT), H5Dclose);
            H5 space(H5Dget_space(data.id), H5Sclose);
            if (H5Sget_simple_extent_ndims(space.id) != static_cast<int>(expected.size()))
                throw std::runtime_error("HDF5 dataset rank mismatch");
            std::vector<hsize_t> dims(expected.size());
            ensure(H5Sget_simple_extent_dims(space.id, dims.data(), nullptr));
            if (dims != expected)
                throw std::runtime_error("HDF5 dataset shape mismatch");
            std::size_t count = 1;
            for (auto n: dims) {
                if (n > 50000000 || count > 50000000 / std::max<hsize_t>(1, n))
                    throw std::runtime_error("HDF5 dataset too large");
                count *= n;
            }
            std::vector<double> value(count);
            if (count)
                ensure(H5Dread(data.id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               value.data()));
            for (double x: value)
                if (!std::isfinite(x))
                    throw std::runtime_error("Nonfinite HDF5 numeric value");
            return value;
        }

        json antenna_json(const Antenna &a) {
            json j = {
                {"kind", a.kind},
                {"orientation", a.orientation},
                {"isotropic", a.isotropic},
                {"shortDipole", a.shortDipole},
                {"dipole", a.dipole},
                {"aperture", a.aperture},
                {"receivePowerReference", a.receivePowerReference},
                {"receiveCalibration", a.receiveCalibration}
            };
            if (a.pattern) {
                const auto &p = *a.pattern;
                j["pattern"] = {
                    {"frequenciesHz", p.frequenciesHz},
                    {"thetaDegrees", p.thetaDegrees},
                    {"phiDegrees", p.phiDegrees},
                    {"radiatedPower", p.radiatedPower},
                    {"acceptedPower", p.acceptedPower},
                    {"stimulatedPower", p.stimulatedPower},
                    {"exportPosition", p.exportPosition},
                    {"exportOrientation", p.exportOrientation},
                    {"provenance", p.provenance},
                    {"frequencyIndependent", p.frequencyIndependent},
                    {"inputConvention", p.inputConvention},
                    {"coefficientScale", p.coefficientScale}
                };
            }
            return j;
        }

        void write_antenna(hid_t file, const std::string &base, const Antenna &a,
                           const Vec3 &position) {
            write_text(file, base + "/definition",
                       {{"position", position}, {"antenna", antenna_json(a)}});
            if (!a.pattern)
                return;
            group(file, base + "/pattern");
            const auto &p = *a.pattern;
            std::vector<double> values;
            values.reserve(p.coefficients.size() * 6);
            for (const auto &v: p.coefficients)
                for (int k = 0; k < 3; ++k) {
                    values.push_back(v[k].real());
                    values.push_back(v[k].imag());
                }
            write_numbers(file, base + "/pattern/coefficients", values,
                          {
                              p.frequenciesHz.size(), p.phiDegrees.size(),
                              p.thetaDegrees.size(), 3, 2
                          });
            write_numbers(file, base + "/pattern/frequency_hz", p.frequenciesHz,
                          {p.frequenciesHz.size()});
            write_numbers(file, base + "/pattern/theta_deg", p.thetaDegrees,
                          {p.thetaDegrees.size()});
            write_numbers(file, base + "/pattern/phi_deg", p.phiDegrees,
                          {p.phiDegrees.size()});
        }

        std::pair<Antenna, Vec3> read_antenna(hid_t file, const std::string &base) {
            const auto document = read_text(file, base + "/definition");
            const auto &j = document.at("antenna");
            Antenna a;
            j.at("kind").get_to(a.kind);
            j.at("orientation").get_to(a.orientation);
            j.at("isotropic").get_to(a.isotropic);
            j.at("shortDipole").get_to(a.shortDipole);
            j.at("dipole").get_to(a.dipole);
            j.at("aperture").get_to(a.aperture);
            j.at("receivePowerReference").get_to(a.receivePowerReference);
            j.at("receiveCalibration").get_to(a.receiveCalibration);
            if (j.contains("pattern")) {
                auto p = std::make_shared<FarfieldData>();
                const auto &v = j.at("pattern");
                v.at("frequenciesHz").get_to(p->frequenciesHz);
                v.at("thetaDegrees").get_to(p->thetaDegrees);
                v.at("phiDegrees").get_to(p->phiDegrees);
                v.at("radiatedPower").get_to(p->radiatedPower);
                v.at("acceptedPower").get_to(p->acceptedPower);
                v.at("stimulatedPower").get_to(p->stimulatedPower);
                v.at("exportPosition").get_to(p->exportPosition);
                v.at("exportOrientation").get_to(p->exportOrientation);
                v.at("provenance").get_to(p->provenance);
                v.at("frequencyIndependent").get_to(p->frequencyIndependent);
                v.at("inputConvention").get_to(p->inputConvention);
                v.at("coefficientScale").get_to(p->coefficientScale);
                const auto values =
                        read_numbers(file, base + "/pattern/coefficients",
                                     {
                                         p->frequenciesHz.size(), p->phiDegrees.size(),
                                         p->thetaDegrees.size(), 3, 2
                                     });
                p->coefficients.resize(values.size() / 6);
                std::size_t index = 0;
                for (auto &field: p->coefficients)
                    for (int k = 0; k < 3; ++k) {
                        field[k] = {values[index], values[index + 1]};
                        index += 2;
                    }
                a.pattern = p;
            }
            return {a, document.at("position").get<Vec3>()};
        }

        json scene_json(const Scene &scene) {
            json j = json::array();
            for (std::size_t i = 0; i < scene.surfaces().size(); ++i) {
                const auto &surface = scene.surfaces()[i];
                std::visit(
                    [&](const auto &shape) {
                        j.push_back({
                            {"type", surface.index()}, {"parameters", shape},
                            {"material", scene.materials()[i]}
                        });
                    },
                    surface);
            }
            return j;
        }

        Scene scene_from_json(const json &j) {
            Scene s;
            for (const auto &entry: j) {
                const auto &p = entry.at("parameters");
                const auto material = entry.at("material").get<Material>();
                switch (entry.at("type").get<int>()) {
                    case 0:
                        s.add(p.get<Rectangle>(), material);
                        break;
                    case 1:
                        s.add(p.get<Disk>(), material);
                        break;
                    case 2:
                        s.add(p.get<Sphere>(), material);
                        break;
                    case 3:
                        s.add(p.get<Cylinder>(), material);
                        break;
                    case 4:
                        s.add(p.get<Triangle>(), material);
                        break;
                    default:
                        throw std::runtime_error("Unknown geometry type");
                }
            }
            return s;
        }

        void parent_directory(const std::filesystem::path &path) {
            if (!path.parent_path().empty())
                std::filesystem::create_directories(path.parent_path());
        }

        std::vector<std::string> split(const std::string &line) {
            std::vector<std::string> out;
            std::istringstream in(line);
            std::string value;
            while (std::getline(in, value, ','))
                out.push_back(value);
            return out;
        }

        double numeric(const std::string &word) {
            std::size_t end = 0;
            const double value = std::stod(word, &end);
            if (end != word.size() || !std::isfinite(value))
                throw std::runtime_error("Invalid impulse CSV number");
            return value;
        }
    } // namespace

    void save_h5(const SimulationResult &result,
                 const std::filesystem::path &path) {
        parent_directory(path);
        H5 file(
            H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT),
            H5Fclose);
        for (const std::string name:
             {
                 "/meta", "/settings", "/scene", "/transmitter", "/receiver", "/paths",
                 "/candidates", "/response", "/diagnostics"
             })
            group(file.id, name);
        write_text(file.id, "/meta/schema",
                   {
                       {"name", "RayTracing"},
                       {"version", schemaVersion},
                       {"solverVersion", RAYTRACING_VERSION},
                       {"lengthUnit", "m"},
                       {"delayUnit", "s"},
                       {"fieldUnit", "V/m"},
                       {"phaseConvention", "exp(-i omega t)"},
                       {"speedOfLight", constants::speedOfLight}
                   });
        write_text(file.id, "/settings/definition", result.settings);
        write_text(file.id, "/scene/geometry", scene_json(result.scene));
        write_antenna(file.id, "/transmitter", result.transmitter.antenna,
                      result.transmitter.position);
        write_antenna(file.id, "/receiver", result.receiver.antenna,
                      result.receiver.position);
        write_text(file.id, "/paths/records", result.rays);
        write_text(file.id, "/candidates/records", result.candidates);
        write_text(file.id, "/response/definition", result.impulseResponse);
        write_text(file.id, "/diagnostics/definition",
                   {
                       {"launchedRays", result.launchedRays},
                       {"traceStatusCounts", result.traceStatusCounts},
                       {"responseRayIndices", result.responseRayIndices}
                   });
        std::vector<double> vertices, offsets{0}, fields, lengths, spreading;
        for (const auto &ray: result.rays) {
            const auto append = [&](const Vec3 &p) {
                for (int k = 0; k < 3; ++k)
                    vertices.push_back(p[k]);
            };
            append(ray.refinement.transmitterPosition_);
            for (const auto &hit: ray.refinement.geometry_.reflections_)
                append(hit.position_);
            if (ray.refinement.geometry_.receiver_)
                append(ray.refinement.geometry_.receiver_->closestPoint_);
            offsets.push_back(vertices.size() / 3);
            lengths.push_back(ray.refinement.geometry_.receiver_
                                  ? ray.refinement.geometry_.receiver_->pathDistance_
                                  : 0);
            for (int k = 0; k < 3; ++k) {
                const Complex value =
                        ray.field ? ray.field->receiverField_[k] : Complex(0);
                fields.push_back(value.real());
                fields.push_back(value.imag());
            }
            for (int r = 0; r < 2; ++r)
                for (int c = 0; c < 2; ++c)
                    spreading.push_back(ray.spreading.jacobian_(r, c));
        }
        write_numbers(file.id, "/paths/vertices", vertices, {vertices.size() / 3, 3});
        write_numbers(file.id, "/paths/vertex_offsets", offsets, {offsets.size()});
        write_numbers(file.id, "/paths/length_m", lengths, {lengths.size()});
        write_numbers(file.id, "/paths/receiver_fields", fields,
                      {result.rays.size(), 3, 2});
        write_numbers(file.id, "/paths/spreading_jacobians", spreading,
                      {result.rays.size(), 2, 2});
        std::vector<double> taps;
        for (const auto &t: result.impulseResponse.taps_) {
            taps.push_back(t.delaySeconds_);
            taps.push_back(t.coefficient_.real());
            taps.push_back(t.coefficient_.imag());
        }
        write_numbers(file.id, "/response/taps", taps,
                      {result.impulseResponse.taps_.size(), 3});
        ensure(H5Fflush(file.id, H5F_SCOPE_GLOBAL));
    }

    SimulationResult load_h5(const std::filesystem::path &path) {
        H5 file(H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
                H5Fclose);
        const auto meta = read_text(file.id, "/meta/schema");
        if (meta.at("name") != "RayTracing" || meta.at("version") != schemaVersion ||
            meta.at("phaseConvention") != "exp(-i omega t)")
            throw std::runtime_error(
                "Unsupported simulation HDF5 schema or phase convention");
        SimulationResult result;
        result.settings =
                read_text(file.id, "/settings/definition").get<SolverConfig>();
        result.scene = scene_from_json(read_text(file.id, "/scene/geometry"));
        auto tx = read_antenna(file.id, "/transmitter"),
                rx = read_antenna(file.id, "/receiver");
        result.transmitter = {tx.second, tx.first};
        result.receiver = {rx.second, rx.first};
        result.transmitter.antenna.validate(result.settings.frequencyHz,
                                            result.settings.medium);
        result.receiver.antenna.validate(result.settings.frequencyHz,
                                         result.settings.medium);
        result.rays =
                read_text(file.id, "/paths/records").get<std::vector<SolvedRay> >();
        result.candidates =
                read_text(file.id, "/candidates/records").get<std::vector<Candidate> >();
        result.impulseResponse =
                read_text(file.id, "/response/definition").get<ImpulseResponse>();
        const auto diag = read_text(file.id, "/diagnostics/definition");
        diag.at("launchedRays").get_to(result.launchedRays);
        diag.at("traceStatusCounts").get_to(result.traceStatusCounts);
        diag.at("responseRayIndices").get_to(result.responseRayIndices);
        if (result.traceStatusCounts.size() != traceStatusCount)
            throw std::runtime_error("Invalid saved trace status counts");
        for (auto index: result.responseRayIndices)
            if (index >= result.rays.size() || !result.rays[index].field)
                throw std::runtime_error("Invalid saved response path index");
        for (const auto &ray: result.rays) {
            if (ray.candidateIndex >= result.candidates.size())
                throw std::runtime_error("Invalid saved candidate index");
            for (const auto &hit: ray.refinement.geometry_.reflections_)
                if (hit.surfaceIndex_ >= result.scene.surfaces().size())
                    throw std::runtime_error("Invalid saved surface index");
        }
        for (const auto &tap: result.impulseResponse.taps_)
            for (auto index: tap.pathIndices_)
                if (index >= result.responseRayIndices.size())
                    throw std::runtime_error("Invalid saved impulse path index");
        return result;
    }

    void save_impulse_csv(const ImpulseResponse &response,
                          const std::filesystem::path &path) {
        parent_directory(path);
        std::ofstream out(path);
        if (!out)
            throw std::runtime_error("Cannot create impulse CSV");
        out << "# RayTracing impulse response 1\n# "
                << json({
                    {"frequencyHz", response.frequencyHz_},
                    {"phaseConvention", "exp(-i omega t)"},
                    {"delayToleranceSeconds", response.delayToleranceSeconds_},
                    {"receiverModel", response.receiverModel_},
                    {"receiverPolarization", response.receiverPolarization_},
                    {
                        "normalization",
                        "common source reference; reciprocal Rx pattern"
                    }
                })
                .dump()
                << '\n';
        out << std::setprecision(17)
                << "delay_s,gain_real,gain_imag,magnitude,phase_deg,path_count,path_"
                "indices,Ex_real,Ex_imag,Ey_real,Ey_imag,Ez_real,Ez_imag\n";
        double peak = 0;
        for (const auto &t: response.taps_)
            peak = std::max(peak, std::abs(t.coefficient_));
        for (const auto &t: response.taps_) {
            const double magnitude = std::abs(t.coefficient_);
            out << t.delaySeconds_ << ',' << t.coefficient_.real() << ','
                    << t.coefficient_.imag() << ',' << magnitude << ',';
            if (magnitude > peak * 1e-12)
                out << std::arg(t.coefficient_) * 180 / constants::pi;
            out << ',' << t.pathIndices_.size() << ',';
            for (std::size_t i = 0; i < t.pathIndices_.size(); ++i) {
                if (i)
                    out << ';';
                out << t.pathIndices_[i];
            }
            for (int k = 0; k < 3; ++k)
                out << ',' << t.normalizedField_[k].real() << ','
                        << t.normalizedField_[k].imag();
            out << '\n';
        }
        out.close();
        if (!out)
            throw std::runtime_error("Failed writing impulse CSV");
    }

    ImpulseResponse load_impulse_csv(const std::filesystem::path &path) {
        std::ifstream input(path);
        if (!input)
            throw std::runtime_error("Cannot open impulse CSV");
        ImpulseResponse response;
        std::string line;
        bool header = false, version = false, metadata = false;
        while (std::getline(input, line)) {
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            if (line.empty())
                continue;
            if (line.starts_with("# RayTracing impulse response")) {
                if (line != "# RayTracing impulse response 1")
                    throw std::runtime_error("Unsupported impulse CSV version");
                version = true;
                continue;
            }
            if (line.starts_with("# {")) {
                const auto meta = json::parse(line.substr(2));
                if (meta.at("phaseConvention") != "exp(-i omega t)")
                    throw std::runtime_error("Unsupported CSV phase convention");
                meta.at("frequencyHz").get_to(response.frequencyHz_);
                meta.at("delayToleranceSeconds").get_to(response.delayToleranceSeconds_);
                meta.at("receiverModel").get_to(response.receiverModel_);
                meta.at("receiverPolarization").get_to(response.receiverPolarization_);
                if (!std::isfinite(response.frequencyHz_) || response.frequencyHz_ < 0 ||
                    !std::isfinite(response.delayToleranceSeconds_) ||
                    response.delayToleranceSeconds_ < 0 ||
                    !response.receiverPolarization_.allFinite())
                    throw std::runtime_error("Invalid impulse CSV metadata");
                metadata = true;
                continue;
            }
            if (line.starts_with('#'))
                continue;
            if (!header) {
                if (line !=
                    "delay_s,gain_real,gain_imag,magnitude,phase_deg,path_count,path_"
                    "indices,Ex_real,Ex_imag,Ey_real,Ey_imag,Ez_real,Ez_imag")
                    throw std::runtime_error("Unrecognized impulse CSV columns");
                header = true;
                continue;
            }
            const auto cols = split(line);
            if (cols.size() != 13)
                throw std::runtime_error("Malformed impulse CSV row");
            ImpulseTap t;
            t.delaySeconds_ = numeric(cols[0]);
            t.coefficient_ = {numeric(cols[1]), numeric(cols[2])};
            if (t.delaySeconds_ < 0 ||
                (!response.taps_.empty() &&
                 t.delaySeconds_ < response.taps_.back().delaySeconds_))
                throw std::runtime_error(
                    "Impulse CSV delays must be nonnegative and sorted");
            std::istringstream ids(cols[6]);
            std::string id;
            while (std::getline(ids, id, ';')) {
                const double value = numeric(id);
                if (value < 0 || value > 1e12 || std::floor(value) != value)
                    throw std::runtime_error("Invalid CSV path index");
                t.pathIndices_.push_back(static_cast<std::size_t>(value));
            }
            if (numeric(cols[5]) != t.pathIndices_.size())
                throw std::runtime_error("CSV path count mismatch");
            for (int k = 0; k < 3; ++k)
                t.normalizedField_[k] = {
                    numeric(cols[7 + 2 * k]),
                    numeric(cols[8 + 2 * k])
                };
            response.taps_.push_back(std::move(t));
        }
        if (!header || !version || !metadata)
            throw std::runtime_error("Missing impulse CSV header, version or metadata");
        return response;
    }
} // namespace rt
