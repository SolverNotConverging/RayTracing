#include "Spreading.hpp"

#include <Eigen/SVD>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <type_traits>

namespace {
using Bundle = Eigen::Matrix<double, 3, 2>;

// On each straight segment det(X+s V) is quadratic. Count its real zeros
// with multiplicity; exclude the source/segment start and final receiver.
unsigned int segment_foci(const Bundle &x, const Bundle &v, const Vec3 &k,
                          double length, bool includeEnd) {
    const Vec3 u = k.unitOrthogonal(), w = k.cross(u);
    Eigen::Matrix2d a, b;
    a.row(0) = u.transpose() * x; a.row(1) = w.transpose() * x;
    b.row(0) = length * u.transpose() * v; b.row(1) = length * w.transpose() * v;
    const long double c0 = a.determinant(), c2 = b.determinant();
    const long double c1 = a(0,0)*b(1,1)+b(0,0)*a(1,1)-a(0,1)*b(1,0)-b(0,1)*a(1,0);
    const long double scale = std::max({std::abs(c0),std::abs(c1),std::abs(c2)});
    const long double eps = 64 * std::numeric_limits<double>::epsilon();
    const auto inside = [&](long double t) {
        return t > eps && (includeEnd ? t <= 1+eps : t < 1-eps);
    };
    if (scale == 0) throw std::runtime_error("Degenerate ray bundle along a segment");
    if (std::abs(c2) <= eps*scale)
        return std::abs(c1) > eps*scale && inside(-c0/c1) ? 1u : 0u;
    const long double discriminant = c1*c1-4*c2*c0;
    const long double tolerance = eps*(c1*c1+std::abs(4*c2*c0));
    if (discriminant < -tolerance) return 0;
    if (std::abs(discriminant) <= tolerance) return inside(-c1/(2*c2)) ? 2u : 0u;
    const long double q = -0.5L*(c1+std::copysign(std::sqrt(discriminant),c1));
    return static_cast<unsigned int>(inside(q/c2)) + static_cast<unsigned int>(inside(c0/q));
}

unsigned int count_caustics(const SequenceEvaluation &path, const Vec3 &launch,
                           const Vec3 &u, const Vec3 &v, const std::vector<Surface> &surfaces) {
    Bundle x = Bundle::Zero(), d;
    d.col(0)=u; d.col(1)=v;
    Vec3 k=launch;
    unsigned int count=0;
    for (const auto &hit : path.reflections_) {
        count += segment_foci(x,d,k,hit.segmentDistance_,true);
        x += hit.segmentDistance_*d;
        const Vec3 n=hit.normal_;
        const double cosine=n.dot(k);
        if (std::abs(cosine)<1e-12) throw std::runtime_error("Grazing variational reflection");
        const Eigen::RowVector2d dt=-(n.transpose()*x)/cosine;
        const Bundle q=x+k*dt;
        Eigen::Matrix3d curvature=Eigen::Matrix3d::Zero();
        std::visit([&](const auto &shape) {
            using T=std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<T,Sphere>)
                curvature=(Eigen::Matrix3d::Identity()-n*n.transpose())/shape.radius_;
            else if constexpr (std::is_same_v<T,Cylinder>) {
                if (std::abs(n.dot(shape.axis_))<0.5)
                    curvature=(Eigen::Matrix3d::Identity()-shape.axis_*shape.axis_.transpose()-n*n.transpose())/shape.radius_;
            }
        },surfaces[hit.surfaceIndex_]);
        const Bundle dn=curvature*q;
        d=(d-2*(n*(n.transpose()*d+k.transpose()*dn)+cosine*dn)).eval();
        k=(k-2*cosine*n).eval();
        x=q-k*dt;
        if (!x.allFinite() || !d.allFinite()) throw std::runtime_error("Variational transport overflow");
    }
    return count+segment_foci(x,d,k,path.receiver_->finalSegmentDistance_,false);
}
}

SpreadingResult calculate_spreading(
    const Vec3 &transmitterPosition, const Vec3 &launchDirection,
    const Vec3 &receiverPosition, const std::vector<Surface> &surfaces,
    const std::vector<std::size_t> &surfaceSequence, const SpreadingOptions &options) {
    const auto positive = [](double x) { return std::isfinite(x) && x > 0.0; };
    if (!positive(options.angularStep_) || options.angularStep_ > 0.1 ||
        options.maxStepHalvings_ < 0 || !positive(options.receiverTolerance_) ||
        !positive(options.relativeDerivativeTolerance_) ||
        !positive(options.absoluteDerivativeTolerance_) ||
        !positive(options.minimumSingularValue_) ||
        !positive(options.minimumSingularValueRatio_) || options.minimumSingularValueRatio_ >= 1.0 ||
        !positive(options.referenceDistance_) || !std::isfinite(options.tMin_) || options.tMin_ < 0.0)
        throw std::invalid_argument("Invalid spreading options");

    const auto central = evaluate_sequence(transmitterPosition, launchDirection,
                                           receiverPosition, surfaces, surfaceSequence, options.tMin_);
    SpreadingResult result;
    if (central.status_ != SequenceStatus::Valid || !central.receiver_ ||
        central.receiver_->missDistance_ > options.receiverTolerance_)
        return result;

    result.launchU_ = launchDirection.unitOrthogonal();
    result.launchV_ = launchDirection.cross(result.launchU_).normalized();
    const Vec3 planeNormal = central.finalDirection_;
    result.receiverU_ = planeNormal.unitOrthogonal();
    result.receiverV_ = planeNormal.cross(result.receiverU_).normalized();

    const auto coordinates = [&](double alpha, double beta) -> std::optional<Eigen::Vector2d> {
        const Vec3 direction = (launchDirection + alpha * result.launchU_ + beta * result.launchV_).normalized();
        const auto trial = evaluate_sequence(transmitterPosition, direction, receiverPosition,
                                             surfaces, surfaceSequence, options.tMin_);
        // The closest-point visibility interval may differ from the fixed-plane
        // interval. Recheck visibility below against the actual plane distance.
        if ((trial.status_ != SequenceStatus::Valid && trial.status_ != SequenceStatus::FinalSegmentBlocked) ||
            trial.reflections_.size() != surfaceSequence.size())
            return std::nullopt;
        const double denominator = trial.finalDirection_.dot(planeNormal);
        if (denominator <= 1e-8) return std::nullopt;
        const double distance = (receiverPosition - trial.finalOrigin_).dot(planeNormal) / denominator;
        if (!std::isfinite(distance) || distance < 0.0) return std::nullopt;
        if (distance > options.tMin_ && nearest_surface(trial.finalOrigin_, trial.finalDirection_,
                                                        surfaces, options.tMin_, distance))
            return std::nullopt;
        const Vec3 offset = trial.finalOrigin_ + distance * trial.finalDirection_ - receiverPosition;
        Eigen::Vector2d value{result.receiverU_.dot(offset), result.receiverV_.dot(offset)};
        if (!value.allFinite()) return std::nullopt;
        return value;
    };
    const auto derivative = [&](double h) -> std::optional<Eigen::Matrix2d> {
        const auto pu = coordinates(h, 0), mu = coordinates(-h, 0);
        const auto pv = coordinates(0, h), mv = coordinates(0, -h);
        if (!pu || !mu || !pv || !mv) return std::nullopt;
        Eigen::Matrix2d j;
        j.col(0) = (*pu - *mu) / (2.0 * h);
        j.col(1) = (*pv - *mv) / (2.0 * h);
        if (!j.allFinite()) return std::nullopt;
        return j;
    };

    double h = options.angularStep_;
    result.status_ = SpreadingStatus::InvalidPerturbation;
    for (int attempt = 0; attempt <= options.maxStepHalvings_; ++attempt, h *= 0.5) {
        const auto coarse = derivative(h), fine = derivative(h * 0.5);
        if (!coarse || !fine) continue;
        result.status_ = SpreadingStatus::UnstableDerivative;
        result.jacobian_ = *fine;
        result.angularStep_ = h * 0.5;
        result.derivativeDifference_ = (*coarse - *fine).norm();
        result.singularValues_ = Eigen::JacobiSVD<Eigen::Matrix2d>(*fine).singularValues();
        result.areaPerSolidAngle_ = std::abs(fine->determinant());
        if (result.derivativeDifference_ > options.absoluteDerivativeTolerance_ +
            options.relativeDerivativeTolerance_ * fine->norm())
            continue;
        // A large singular value must not hide uncertainty in a small one.
        const Eigen::Vector2d coarseSingularValues = Eigen::JacobiSVD<Eigen::Matrix2d>(*coarse).singularValues();
        bool stable = true;
        for (int i = 0; i < 2; ++i)
            stable = stable && std::abs(coarseSingularValues[i] - result.singularValues_[i]) <=
                     options.absoluteDerivativeTolerance_ + options.relativeDerivativeTolerance_ * result.
                     singularValues_[i];
        if (!stable) continue;
        if (result.singularValues_[1] <= options.minimumSingularValue_ ||
            result.singularValues_[1] <= options.minimumSingularValueRatio_ * result.singularValues_[0]) {
            result.status_ = SpreadingStatus::Caustic;
            return result;
        }
        const double factor = options.referenceDistance_ / std::sqrt(result.areaPerSolidAngle_);
        if (!std::isfinite(result.areaPerSolidAngle_) || !positive(factor)) return result;
        result.fieldFactor_ = factor;
        result.causticCount_ = count_caustics(central, launchDirection, result.launchU_, result.launchV_, surfaces);
        result.status_ = SpreadingStatus::Valid;
        return result;
    }
    return result;
}

const char *spreading_status_name(SpreadingStatus status) {
    switch (status) {
        case SpreadingStatus::Valid: return "valid";
        case SpreadingStatus::InvalidCentralPath: return "invalid central path";
        case SpreadingStatus::InvalidPerturbation: return "invalid neighbouring ray";
        case SpreadingStatus::UnstableDerivative: return "unstable derivative";
        case SpreadingStatus::Caustic: return "near receiver caustic";
    }
    return "unknown";
}
