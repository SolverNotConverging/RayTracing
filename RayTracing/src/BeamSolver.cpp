#include <RayTracing/BeamSolver.hpp>
#include <RayTracing/Constants.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace rt::beam {
namespace {
constexpr double pi = rt::constants::pi;
constexpr double c0 = rt::constants::speedOfLight;
constexpr double eta0 = rt::constants::freeSpaceImpedance;
bool finite(Complex c) { return std::isfinite(c.real()) && std::isfinite(c.imag()); }
void require(bool condition, const char *message) { if (!condition) throw std::invalid_argument(message); }
void frequency_valid(double f) { require(std::isfinite(f) && f > 0, "Frequency must be positive and finite"); }
double sign(double x) { return (x > 0) - (x < 0); }
}
void Material::validate() const {
    require(kind == "dielectric" || kind == "metal" || kind == "pec", "Unknown material kind");
    require(finite(epsilon) && epsilon.real() > 0 && epsilon.imag() >= 0 && std::isfinite(mu) && mu > 0 &&
        std::isfinite(conductivity) && conductivity >= 0, "Require passive finite material parameters");
}
Complex Material::permittivity(double f) const {
    validate(); frequency_valid(f);
    require(kind != "pec", "PEC has no finite permittivity");
    return epsilon + Complex(0, conductivity / (2 * pi * f * rt::constants::vacuumPermittivity));
}
Complex Material::index(double f) const { return std::sqrt(permittivity(f) * mu); }
std::array<Complex, 4> fresnel(const Material &incident, const Material &outgoing, double c, double f) {
    incident.validate(); outgoing.validate(); frequency_valid(f);
    require(c >= 0 && c <= 1 && incident.kind == "dielectric", "Require dielectric incident medium and cosine in [0,1]");
    const auto n1 = incident.index(f);
    if (outgoing.kind == "pec") return {-1., 1., 0., 0.};
    const auto e1 = incident.permittivity(f), e2 = outgoing.permittivity(f), n2 = outgoing.index(f);
    const auto q1 = n1 * c;
    auto q2 = std::sqrt(n2*n2 - n1*n1*(1-c*c));
    if (q2.imag() < 0 || (q2.imag() == 0 && q2.real() < 0)) q2 = -q2;
    auto ratio = [](Complex a, Complex b) -> Complex { return std::abs(a+b) > 0 ? (a-b)/(a+b) : Complex(0); };
    const auto rs = ratio(outgoing.mu*q1, incident.mu*q2), rp = ratio(e2*q1, e1*q2);
    if (outgoing.kind == "metal") return {rs, rp, 0., 0.};
    return {rs, rp, 1.+rs, n1/n2*(1.+rp)};
}
Vec3 unit(const Vec3 &v) {
    require(v.allFinite() && v.norm() > 0, "Expected a finite nonzero direction");
    return v.normalized();
}
Frame frame(const Vec3 &direction) {
    const auto d = unit(direction);
    Eigen::Index i; d.cwiseAbs().minCoeff(&i);
    const Vec3 u = Vec3::Unit(i).cross(d).normalized();
    Frame result; result.col(0) = u; result.col(1) = d.cross(u); return result;
}
Vec3C polarization(const Vec3C &p, const Vec3 &d) {
    require(p.allFinite() && p.norm() > 0, "Polarization must be a finite nonzero complex 3-vector");
    require(std::abs((d.transpose().cast<Complex>() * p)(0)) <= 1e-10 * p.norm(), "Polarization must be transverse");
    return p / p.norm();
}
void Box::validate() const {
    require(center.allFinite() && halfLengths.allFinite() && (halfLengths.array() > 0).all() && axes.allFinite(), "Box half lengths must be positive and finite");
    require((axes.transpose()*axes - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() <= 1e-12 && axes.determinant() > 0,
        "Box axes must be a right-handed orthonormal frame");
}
std::vector<Surface> faces(const Shape &shape) {
    return std::visit([](const auto &s) -> std::vector<Surface> {
        if constexpr (std::is_same_v<std::decay_t<decltype(s)>, Box>) {
            s.validate(); std::vector<Surface> result;
            for (int i=0;i<3;++i) {
                const int j=(i+1)%3, k=(i+2)%3;
                for (double sign : {-1.,1.}) result.push_back(Rectangle{s.center+sign*s.halfLengths[i]*s.axes.col(i),
                    s.axes.col(j), sign*s.axes.col(k), s.halfLengths[j], s.halfLengths[k]});
            }
            return result;
        } else return {s};
    }, shape);
}
bool contains(const Shape &shape, const Vec3 &p) {
    require(p.allFinite(), "Expected a finite 3-vector");
    return std::visit([&](const auto &s) -> bool {
        using T = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<T,Box>) return ((s.axes.transpose()*(p-s.center)).cwiseAbs().array() < s.halfLengths.array()).all();
        else if constexpr (std::is_same_v<T,Sphere>) return (p-s.center_).norm() < s.radius_;
        else if constexpr (std::is_same_v<T,Cylinder>) {
            const Vec3 v=p-s.center_; const double a=v.dot(s.axis_);
            return std::abs(a) < s.halfLength_ && (v-a*s.axis_).norm() < s.radius_;
        } else return false;
    }, shape);
}
double edge_distance(const Surface &shape, const Vec3 &p) {
    return std::visit([&](const auto &s) -> double {
        using T = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<T,Rectangle>) {
            const Vec3 v=p-s.center_; return std::min(s.halfWidth_-std::abs(v.dot(s.u_)), s.halfHeight_-std::abs(v.dot(s.v_)));
        } else if constexpr (std::is_same_v<T,Disk>) return s.radius_-(p-s.center_).norm();
        else if constexpr (std::is_same_v<T,Cylinder>) {
            const Vec3 v=p-s.center_; const double axial=v.dot(s.axis_);
            return std::hypot(s.halfLength_-std::abs(axial), s.radius_-(v-axial*s.axis_).norm());
        } else if constexpr (std::is_same_v<T,Triangle>) {
            return std::min({(s.b_-s.a_).cross(p-s.a_).norm()/(s.b_-s.a_).norm(),
                (s.c_-s.b_).cross(p-s.b_).norm()/(s.c_-s.b_).norm(), (s.a_-s.c_).cross(p-s.c_).norm()/(s.a_-s.c_).norm()});
        } else return std::numeric_limits<double>::infinity();
    }, shape);
}
Eigen::Matrix3d curvature(const Surface &shape, const Vec3 &p, const Vec3 &normal) {
    return std::visit([&](const auto &s) -> Eigen::Matrix3d {
        using T = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<T,Sphere>) {
            const Vec3 outward=(p-s.center_)/s.radius_;
            return sign(outward.dot(normal))*(Eigen::Matrix3d::Identity()-outward*outward.transpose())/s.radius_;
        } else if constexpr (std::is_same_v<T,Cylinder>) {
            const Vec3 v=p-s.center_, radial=v-v.dot(s.axis_)*s.axis_;
            if (radial.norm()>0 && std::abs(normal.dot(s.axis_))<.5) {
                const Vec3 outward=radial.normalized();
                return sign(outward.dot(normal))*(Eigen::Matrix3d::Identity()-s.axis_*s.axis_.transpose()-outward*outward.transpose())/s.radius_;
            }
        }
        return Eigen::Matrix3d::Zero();
    }, shape);
}
Scene::Scene(Material background, double tolerance) : background_(std::move(background)), tolerance_(tolerance) {
    background_.validate(); require(background_.kind == "dielectric", "Background must be a dielectric");
    require(std::isfinite(tolerance) && tolerance > 0, "Tolerance must be positive");
}
std::size_t Scene::add(Shape shape, std::optional<Material> material, bool volume, std::string terminate) {
    if (volume) require(std::holds_alternative<Box>(shape) || std::holds_alternative<Sphere>(shape) ||
        (std::holds_alternative<Cylinder>(shape) && std::get<Cylinder>(shape).capped_), "Volumes must be boxes, spheres, or capped cylinders");
    require(terminate.empty() || terminate == "entry" || (volume && terminate == "exit"), "Use entry/exit for volumes; surfaces terminate on any hit");
    if (terminate.empty()) {
        require(material.has_value(), "Expected Material"); material->validate();
        require(volume || material->kind != "dielectric", "Dielectrics require a volume; surfaces must be PEC or metal");
    }
    const auto fs=faces(shape); for (const auto &f : fs) validate_surface(f);
    const auto id=objects_.size();
    objects_.push_back({std::move(shape),std::move(material),volume,std::move(terminate)});
    for (const auto &f : fs) faces_.emplace_back(id,f);
    return id;
}
void Scene::validate() const { for (const auto &[id,f] : faces_) validate_surface(f); }
int Scene::region_at(const Vec3 &p) const {
    require(p.allFinite(), "Expected a finite 3-vector");
    for (std::size_t i=objects_.size();i>0;--i) {
        const auto &o=objects_[i-1]; if (o.volume && o.terminate.empty() && contains(o.shape,p)) return static_cast<int>(i-1);
    }
    return -1;
}
const Material &Scene::material(int region) const {
    if (region == -1) return background_;
    require(region >= 0 && static_cast<std::size_t>(region) < objects_.size() && objects_[region].material.has_value(), "Invalid material region");
    return *objects_[region].material;
}
bool Scene::starts_terminated(const Vec3 &p) const {
    require(p.allFinite(), "Expected a finite 3-vector");
    return std::any_of(objects_.begin(),objects_.end(),[&](const auto &o) { return o.volume && o.terminate == "entry" && contains(o.shape,p); });
}
std::optional<Event> Scene::next_event(const Vec3 &origin, const Vec3 &direction, double maxDistance) const {
    require(origin.allFinite() && direction.allFinite() && std::abs(direction.norm()-1)<1e-10 && std::isfinite(maxDistance) && maxDistance >= 0,
        "Require finite origin, unit direction and nonnegative distance");
    struct Hit { double distance; std::size_t id; Surface shape; Vec3 normal; };
    Vec3 current=origin; double travelled=0;
    // Iterate invisible crossings; no recursion limit for deeply nested scenes.
    while (maxDistance-travelled > tolerance_) {
        std::vector<Hit> hits;
        for (const auto &[id,f] : faces_) {
            const auto hit=std::visit([&](const auto &s) { return intersect(current,direction,s,tolerance_,maxDistance-travelled); }, f);
            if (hit) hits.push_back({hit->distance_,id,f,hit->normal_});
        }
        if (hits.empty()) return std::nullopt;
        std::stable_sort(hits.begin(),hits.end(),[](const auto &a,const auto &b) { return a.distance < b.distance; });
        const double distance=hits[0].distance, probe=2*tolerance_;
        hits.erase(std::remove_if(hits.begin(),hits.end(),[&](const auto &h) { return h.distance-distance>tolerance_; }),hits.end());
        const Vec3 p=current+distance*direction;
        const int before=region_at(p-probe*direction), after=region_at(p+probe*direction);
        auto event=[&](const Hit &h, const char *kind) -> Event { return {travelled+distance,p,h.normal,h.shape,h.id,before,after,kind}; };
        for (const auto &h : hits) {
            const auto &o=objects_[h.id];
            if (!o.terminate.empty()) {
                const bool was=contains(o.shape,p-probe*direction), now=contains(o.shape,p+probe*direction);
                if (!o.volume || (o.terminate=="entry" && !was && now) || (o.terminate=="exit" && was && !now)) return event(h,"terminate");
            }
        }
        std::vector<Hit> active;
        for (const auto &h : hits) {
            const auto &o=objects_[h.id];
            if (o.terminate.empty() && (!o.volume || static_cast<int>(h.id)==before || static_cast<int>(h.id)==after)) active.push_back(h);
        }
        if (!active.empty()) {
            const Hit *winner=&active[0]; bool sheet=false;
            for (const auto &h : active) {
                const bool isSheet=!objects_[h.id].volume;
                if ((isSheet && !sheet) || (isSheet==sheet && h.id>winner->id)) winner=&h;
                sheet=sheet || isSheet;
            }
            for (const auto &h : active)
                if (std::abs(winner->normal.dot(h.normal))<1-1e-10 || edge_distance(h.shape,p)<=tolerance_) return event(*winner,"ambiguous");
            if (sheet || !(material(before)==material(after))) return event(*winner,"interface");
        }
        current=p+probe*direction; travelled+=distance+probe;
    }
    return std::nullopt;
}
void State::validate() const {
    require(position.allFinite() && direction.allFinite() && std::abs(direction.norm()-1)<1e-10 && basis.allFinite() &&
        (basis.transpose()*basis-Eigen::Matrix2d::Identity()).norm()<1e-10 && (basis.transpose()*direction).norm()<1e-10,
        "Invalid beam position, direction or transverse frame");
    require(curvature.allFinite() && field.allFinite() && std::isfinite(power) && power>=0 && std::isfinite(distance) && distance>=0 &&
        std::isfinite(opticalPath), "Invalid beam curvature, field or path state");
    require((curvature-curvature.transpose()).norm()<=1e-10*std::max(1.,curvature.norm()), "Beam curvature must be symmetric");
}
State propagate(const State &state, double distance, const Material &material, double frequency) {
    state.validate();
    require(std::isfinite(distance) && distance>=0, "Propagation distance must be nonnegative and finite");
    const Complex n=material.index(frequency);
    const double k0=2*pi*frequency/c0;
    const auto &h=state.curvature;
    const Complex halfTrace=h.trace()/2.;
    const Complex delta=std::sqrt((h(0,0)-h(1,1))*(h(0,0)-h(1,1))/4.+h(0,1)*h(1,0));
    const Complex a=1.+distance*(halfTrace+delta), b=1.+distance*(halfTrace-delta);
    require(std::min(std::abs(a),std::abs(b))>=1e-14, "Singular geometrical plane-wave focus; use a finite GaussianBeam");
    const Complex denominator=std::sqrt(a)*std::sqrt(b);
    State result=state;
    result.position+=distance*state.direction;
    result.curvature=h*(Mat2C::Identity()+distance*h).inverse();
    result.field=state.field/denominator*std::exp(Complex(0,1)*k0*n*distance);
    result.power*=std::exp(-2*k0*n.imag()*distance);
    result.distance+=distance; result.opticalPath+=n.real()*distance;
    return result;
}
Vec3C Segment::field_at(const Vec3 &point) const {
    require(point.allFinite(), "Expected a finite 3-vector");
    const Vec3 offset=point-state.position;
    const double s=offset.dot(state.direction);
    if (s<0 || s>=length) return Vec3C::Zero();
    const auto at=propagate(state,s,material,frequency);
    const Eigen::Vector2cd x=(state.basis.transpose()*offset).cast<Complex>();
    const double k=2*pi*frequency/c0*material.index(frequency).real();
    return at.field*std::exp(Complex(0,.5)*k*(x.transpose()*at.curvature*x)(0));
}
Eigen::Vector2d Segment::radii(double distance) const {
    const auto at=propagate(state,distance,material,frequency);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(at.curvature.imag());
    const auto values=eig.eigenvalues();
    const double k=2*pi*frequency/c0*material.index(frequency).real();
    if (values.minCoeff()<=0) return Eigen::Vector2d::Constant(std::numeric_limits<double>::infinity());
    return (2/(k*values.array())).sqrt().matrix();
}
State GaussianBeam::launch(const Scene &scene, double frequency) const {
    frequency_valid(frequency);
    const auto d=unit(direction);
    require(waistPosition.allFinite(), "Expected a finite 3-vector");
    const Vec3 origin=waistPosition+launchDistance*d;
    Frame basis=frame(d);
    if (transverseAxis) {
        const auto u=unit(*transverseAxis);
        require(std::abs(u.dot(d))<=1e-10, "Transverse axis must be perpendicular to beam axis");
        basis.col(0)=u; basis.col(1)=d.cross(u);
    }
    require(waistRadius.allFinite() && (waistRadius.array()>0).all(), "Waist radii must be positive");
    require(std::isfinite(power) && power>0 && std::isfinite(phase) && std::isfinite(launchDistance), "Require positive finite power and finite launch distance/phase");
    const int region=scene.region_at(origin); const auto &medium=scene.material(region);
    require(medium.kind=="dielectric", "Cannot launch inside opaque material");
    const auto n=medium.index(frequency);
    require(n.imag()==0, "Launch Gaussian sources in a lossless region; absorption begins at subsequent interfaces");
    const double k=2*pi*frequency/c0*n.real();
    require((2/(k*waistRadius.array())).maxCoeff()<=.2, "Gaussian divergence exceeds 0.2 rad; outside this paraxial model");
    const Eigen::Vector2d zr=k*waistRadius.array().square()/2;
    Mat2C h=Mat2C::Zero(); Complex gouy=1;
    for (int i=0;i<2;++i) {
        const Complex q(launchDistance,-zr[i]); h(i,i)=1./q; gouy*=std::sqrt(Complex(0,-zr[i])/q);
    }
    const double admittance=(n/medium.mu).real()/eta0;
    const double amplitude=std::sqrt(4*power/(pi*waistRadius.prod()*admittance));
    State result; result.position=origin; result.direction=d; result.basis=basis; result.curvature=h;
    result.field=amplitude*gouy*std::exp(Complex(0,k*launchDistance+phase))*beam::polarization(polarization,d);
    result.region=region; result.power=power; return result;
}
std::vector<GaussianBeam> GaussianBeam::beamlets(int samples, double ratio, double extent) const {
    require(samples>=3 && ratio>0 && ratio<1 && std::isfinite(extent) && extent>0, "Require samples>=3, 0<radius_ratio<1 and positive extent");
    const auto d=unit(direction); Frame basis=frame(d);
    if (transverseAxis) {
        const auto u=unit(*transverseAxis);
        require(std::abs(u.dot(d))<=1e-10,"Transverse axis must be perpendicular to beam axis");
        basis.col(0)=u; basis.col(1)=d.cross(u);
    }
    require(waistRadius.allFinite() && (waistRadius.array()>0).all() && waistPosition.allFinite(), "Invalid beam waist");
    const Eigen::Vector2d small=waistRadius*ratio;
    const Eigen::Vector2d spread=(waistRadius.array().square()-small.array().square()).sqrt();
    const Eigen::Vector2d step=2*extent*spread/(samples-1);
    std::vector<GaussianBeam> result;
    for (int i=0;i<samples;++i) for (int j=0;j<samples;++j) {
        const Eigen::Vector2d xy(-extent*spread[0]+i*step[0],-extent*spread[1]+j*step[1]);
        const double weight=step.prod()*(i==0 || i==samples-1 ? .5 : 1)*(j==0 || j==samples-1 ? .5 : 1);
        const double amplitude=weight*waistRadius.prod()/(pi*small.prod()*spread.prod())*std::exp(-(xy.array()/spread.array()).square().sum());
        auto beam=*this; beam.waistPosition+=basis*xy; beam.direction=d; beam.waistRadius=small;
        beam.power=power*amplitude*amplitude*small.prod()/waistRadius.prod(); beam.transverseAxis=basis.col(0);
        result.push_back(std::move(beam));
    }
    return result;
}
State PlaneWave::launch(const Scene &scene, double frequency) const {
    frequency_valid(frequency); const auto d=unit(direction);
    const auto region=scene.region_at(position); const auto &mat=scene.material(region);
    require(mat.kind=="dielectric", "Cannot launch inside opaque material");
    require(finite(amplitude), "Amplitude must be finite");
    State result; result.position=position; result.direction=d; result.basis=frame(d);
    result.field=amplitude*beam::polarization(polarization,d); result.region=region;
    result.power=std::norm(amplitude)*(mat.index(frequency)/mat.mu).real()/(2*eta0); return result;
}
std::vector<State> ApertureSource::launch(const Scene &scene, double frequency) const {
    frequency_valid(frequency);
    require(!frequencyHz || frequency==*frequencyHz,"Sampled antenna source frequency must match solver frequency");
    require(samples>=3 && samples<=41,"Aperture samples per axis must lie in [3,41]");
    require(halfWidths.allFinite() && (halfWidths.array()>0).all() && std::isfinite(beamletRadius) && beamletRadius>0,
        "Aperture widths and beamlet radius must be positive");
    require(std::isfinite(rcond) && rcond>0 && rcond<1,"rcond must be in (0,1)");
    require(static_cast<bool>(field),"Aperture field callback is required");
    const auto d=unit(direction); const auto basis=frame(d); require(position.allFinite(),"Expected a finite 3-vector");
    const int count=samples*samples;
    Eigen::MatrixXd xy(count,2), valuesReal(count,3), valuesImag(count,3);
    std::vector<Vec3> positions; positions.reserve(count);
    double norm=0, longitudinal=0;
    for (int i=0;i<samples;++i) for (int j=0;j<samples;++j) {
        const int row=i*samples+j;
        xy.row(row)=Eigen::Vector2d(-halfWidths[0]+2*halfWidths[0]*i/(samples-1),-halfWidths[1]+2*halfWidths[1]*j/(samples-1));
        const Vec3 p=position+basis*xy.row(row).transpose(); positions.push_back(p);
        const auto e=field(p); require(e.allFinite(),"Aperture field must return finite complex Cartesian vectors");
        valuesReal.row(row)=e.real(); valuesImag.row(row)=e.imag(); norm+=e.squaredNorm();
        longitudinal=std::max(longitudinal,std::abs((d.transpose().cast<Complex>()*e)(0)));
    }
    require(longitudinal<=1e-8*std::max(1e-300,std::sqrt(norm)),"Aperture field must be transverse to the launch axis");
    Eigen::MatrixXd kernel(count,count);
    for (int i=0;i<count;++i) for (int j=0;j<count;++j) kernel(i,j)=std::exp(-(xy.row(i)-xy.row(j)).squaredNorm()/(beamletRadius*beamletRadius));
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(kernel,Eigen::ComputeThinU | Eigen::ComputeThinV); svd.setThreshold(rcond);
    const Eigen::MatrixXd real=svd.solve(valuesReal), imag=svd.solve(valuesImag);
    const auto &medium=scene.material(scene.region_at(position));
    require(medium.kind=="dielectric","Cannot launch inside opaque material");
    require(!vacuumOnly || (medium.epsilon==Complex(1) && medium.mu==1 && medium.conductivity==0),"Antenna aperture adapter requires a vacuum launch medium");
    const auto n=medium.index(frequency); require(n.imag()==0,"Aperture sources require a homogeneous lossless launch plane");
    std::vector<State> states;
    for (int i=0;i<count;++i) {
        require(scene.material(scene.region_at(positions[i]))==medium,"Aperture launch plane crosses a material boundary");
        Vec3C e=real.row(i).transpose().cast<Complex>()+Complex(0,1)*imag.row(i).transpose().cast<Complex>();
        const double amplitude=e.norm(); if (amplitude<=1e-30) continue;
        GaussianBeam beam; beam.waistPosition=positions[i]; beam.direction=d; beam.waistRadius.setConstant(beamletRadius);
        beam.power=pi*beamletRadius*beamletRadius*(n/medium.mu).real()/(4*eta0)*amplitude*amplitude;
        beam.polarization=e/amplitude; beam.transverseAxis=basis.col(0); states.push_back(beam.launch(scene,frequency));
    }
    require(!states.empty(),"Aperture field is zero"); return states;
}
ApertureSource ApertureSource::from_antenna(const rt::Antenna &antenna, const Vec3 &antennaPosition,
    const Vec3 &position, const Vec3 &direction, const Eigen::Vector2d &halfWidths, int samples, double radius, double frequency) {
    frequency_valid(frequency); require(antennaPosition.allFinite(),"Expected a finite 3-vector"); const auto d=unit(direction);
    ApertureSource source; source.position=position; source.direction=d; source.halfWidths=halfWidths;
    source.samples=samples; source.beamletRadius=radius; source.frequencyHz=frequency; source.vacuumOnly=true;
    source.field=[antenna,antennaPosition,d,frequency](const Vec3 &p) -> Vec3C {
        const Vec3 offset=p-antennaPosition; const double distance=offset.norm();
        require(distance>0,"Antenna cannot lie on its sampled aperture"); const Vec3 look=offset/distance;
        require(look.dot(d)>=std::cos(.2),"Antenna aperture exceeds the 0.2 rad paraxial cone");
        const Vec3C e=antenna.farfield(look,frequency)/distance;
        return (e-d.cast<Complex>()*(d.transpose().cast<Complex>()*e)(0))*std::exp(Complex(0,2*pi*frequency/c0*distance));
    };
    return source;
}
std::pair<std::vector<State>,std::string> interact(const State &state, const Event &event, const Scene &scene, double frequency) {
    const auto &d=state.direction; const auto &u=state.basis;
    Vec3 normal=event.normal; if (normal.dot(d)>0) normal=-normal;
    const double c=-normal.dot(d); if (c<1e-10) return {{},"grazing"};
    const auto &incoming=scene.material(state.region);
    const auto &object=scene.objects().at(event.objectId);
    const auto &outgoing=object.volume ? scene.material(event.after) : *object.material;
    const double n1=incoming.index(frequency).real(), n2=outgoing.kind=="pec" ? 1 : outgoing.index(frequency).real();
    const auto coefficients=fresnel(incoming,outgoing,c,frequency);
    Vec3 s=d.cross(normal); s=s.norm()>1e-10 ? s.normalized().eval() : u.col(0).eval();
    const Vec3 p=s.cross(d);
    const Complex es=(s.transpose().cast<Complex>()*state.field)(0), ep=(p.transpose().cast<Complex>()*state.field)(0);
    const Eigen::Matrix3d dnDx=curvature(event.shape,event.position,normal);
    const Eigen::Matrix3d projection=Eigen::Matrix3d::Identity()-d*normal.transpose()/normal.dot(d);
    const Frame dx=projection*u, dn=dnDx*dx;
    std::vector<State> children;
    for (const std::string kind : {"reflect","transmit"}) {
        if (kind=="transmit" && outgoing.kind!="dielectric") continue;
        Vec3 direction; Frame derivativeX, derivativeD; Complex cs,cp; int region;
        if (kind=="reflect") {
            direction=d+2*c*normal; cs=coefficients[0]; cp=coefficients[1];
            derivativeX=-2*(normal*(d.transpose()*dn)+d.dot(normal)*dn);
            derivativeD=(Eigen::Matrix3d::Identity()-2*normal*normal.transpose())*u; region=state.region;
        } else {
            const double eta=n1/n2, ct2=1-eta*eta*(1-c*c);
            if (ct2<=1e-14) {
                if (outgoing.index(frequency).imag()!=0 || incoming.index(frequency).imag()!=0) return {{},"unsupported_lossy_critical"};
                continue;
            }
            const double ct=std::sqrt(ct2), derivative=eta-eta*eta*c/ct;
            direction=eta*d+(eta*c-ct)*normal; cs=coefficients[2]; cp=coefficients[3];
            derivativeX=(eta*c-ct)*dn-normal*(derivative*(d.transpose()*dn));
            derivativeD=eta*u-normal*(derivative*(normal.transpose()*u)); region=event.after;
        }
        direction=unit(direction); const Frame v=frame(direction);
        const Eigen::Matrix2d a=v.transpose()*dx;
        Mat2C h=((v.transpose()*derivativeX).cast<Complex>()+(v.transpose()*derivativeD).cast<Complex>()*state.curvature)*a.inverse().cast<Complex>();
        if ((h-h.transpose()).norm()>1e-7*std::max(1.,h.norm())) throw std::runtime_error("Interface curvature lost reciprocity");
        h=((h+h.transpose())/2.).eval();
        const Vec3C field=cs*es*s.cast<Complex>()+cp*ep*s.cross(direction).cast<Complex>();
        double ratio=field.squaredNorm()/std::max(state.field.squaredNorm(),1e-300);
        if (kind=="transmit") ratio*=(outgoing.index(frequency)/outgoing.mu).real()*std::abs(direction.dot(normal))/((incoming.index(frequency)/incoming.mu).real()*c);
        State child=state; child.position=event.position; child.direction=direction; child.basis=v; child.curvature=h;
        child.field=field; child.region=region; child.power=state.power*ratio;
        child.branch=0; child.parent=-1; child.history.emplace_back(event.objectId,kind,state.region,region);
        children.push_back(std::move(child));
    }
    return {std::move(children),""};
}
void SolverConfig::validate() const {
    frequency_valid(frequencyHz);
    require(std::isfinite(maxDistance) && maxDistance>0,"Frequency and max distance must be positive");
    require(maxInteractions>=0 && maxBranches>=1,"Invalid interaction/branch budget");
    require(minPowerFraction>=0 && minPowerFraction<1 && maxLossRatio>=0 && maxLossRatio<=.01,
        "Require 0 <= power threshold < 1 and 0 <= weak-loss ratio <= 0.01");
}
Vec3C Simulation::field_at(const Vec3 &p) const {
    const auto &mat=scene.material(scene.region_at(p));
    if (mat.kind!="dielectric" || scene.starts_terminated(p)) return Vec3C::Zero();
    Vec3C result=Vec3C::Zero();
    for (const auto &segment : segments) if (segment.material==mat) result+=segment.field_at(p);
    return result;
}
std::vector<Vec3C> Simulation::fields_at(const std::vector<Vec3> &points) const {
    std::vector<Vec3C> values; values.reserve(points.size());
    for (const auto &p : points) values.push_back(field_at(p)); return values;
}
std::map<std::string,double> Simulation::power_balance() const {
    std::map<std::string,double> result{{"launched",launchedPower},{"absorbed",absorbedPower},{"interface_flux_residual",interfaceFluxResidual}};
    for (const auto &t : terminals) result[t.reason]+=t.power; return result;
}
Simulation solve_states(const Scene &scene, const std::vector<State> &states, const SolverConfig &config) {
    config.validate(); scene.validate();
    auto checkLoss=[&](const Material &m) {
        if (m.kind=="dielectric") { const auto n=m.index(config.frequencyHz);
            require(n.imag()/n.real()<=config.maxLossRatio,"Dielectric loss exceeds configured weak-loss limit; strongly absorbing refraction is unsupported"); }
    };
    checkLoss(scene.background()); for (const auto &o : scene.objects()) if (o.material) checkLoss(*o.material);
    require(!states.empty(),"At least one source is required"); require(states.size()<=config.maxBranches,"Source count exceeds branch budget");
    Simulation result; result.config=config; result.scene=scene;
    std::vector<State> queue=states;
    for (std::size_t i=0;i<queue.size();++i) {
        queue[i].validate(); require(scene.material(queue[i].region).kind=="dielectric","Cannot launch inside opaque material");
        queue[i].branch=i; result.launchedPower+=queue[i].power;
    }
    require(std::isfinite(result.launchedPower),"Total launched power must be finite");
    auto warn=[&](const std::string &message) {
        if (std::find(result.warnings.begin(),result.warnings.end(),message)==result.warnings.end()) result.warnings.push_back(message);
    };
    if (queue.size()>1) warn("Coherent source superposition: power_balance sums constituent branch powers; use a field flux integral for physical total power");
    std::size_t allocated=queue.size(), nextId=allocated;
    auto finish=[&](const State &s,const std::string &reason) { result.terminals.push_back({s.branch,reason,s.position,s.power,s.history}); };
    while (!queue.empty()) {
        State state=std::move(queue.back()); queue.pop_back();
        if (scene.starts_terminated(state.position)) { finish(state,"terminated"); continue; }
        if (state.power==0 || state.power<config.minPowerFraction*result.launchedPower) { finish(state,"power_cutoff"); continue; }
        const double remaining=config.maxDistance-state.distance;
        if (remaining<=scene.tolerance()) { finish(state,"distance_limit"); continue; }
        const auto event=scene.next_event(state.position,state.direction,remaining);
        const double length=event ? event->distance : remaining;
        const auto &material=scene.material(state.region);
        Segment segment{state,length,material,config.frequencyHz}; result.segments.push_back(segment);
        State end=propagate(state,length,material,config.frequencyHz); result.absorbedPower+=state.power-end.power;
        if (!event) { finish(end,"distance_limit"); continue; }
        if (event->kind=="ambiguous" || event->kind=="terminate") { finish(end,event->kind=="ambiguous" ? "ambiguous" : "terminated"); continue; }
        if (state.history.size()>=static_cast<std::size_t>(config.maxInteractions)) { finish(end,"interaction_limit"); continue; }
        end.region=event->before;
        const auto widths=segment.radii(length);
        if (widths.allFinite()) {
            const double footprint=3*widths.maxCoeff()/std::max(std::abs(state.direction.dot(event->normal)),1e-10);
            if (edge_distance(event->shape,event->position)<footprint)
                warn("Branch "+std::to_string(state.branch)+": beam footprint approaches object "+std::to_string(event->objectId)+" edge; clipping is unresolved");
        }
        auto [children,reason]=interact(end,*event,scene,config.frequencyHz);
        if (!reason.empty()) { finish(end,reason); continue; }
        double outgoingPower=0; for (const auto &child : children) outgoingPower+=child.power;
        const double residual=end.power-outgoingPower;
        const auto &target=scene.objects()[event->objectId];
        const auto &boundary=target.volume ? scene.material(event->after) : *target.material;
        if (boundary.kind=="metal" || boundary.kind=="pec") result.absorbedPower+=residual;
        else {
            result.interfaceFluxResidual+=residual;
            if (std::abs(residual)>1e-7*std::max(end.power,1e-300)) warn("Lossy-interface branch budgets omit interference flux; inspect interface_flux_residual");
        }
        for (auto &child : children) {
            child.parent=static_cast<std::int64_t>(state.branch); child.branch=nextId++;
            if (allocated>=config.maxBranches) finish(child,"branch_limit");
            else { ++allocated; queue.push_back(std::move(child)); }
        }
    }
    return result;
}
Simulation solve(const Scene &scene, const std::vector<Source> &sources, const SolverConfig &config) {
    config.validate(); std::vector<State> states;
    for (const auto &source : sources) std::visit([&](const auto &s) {
        if constexpr (std::is_same_v<std::decay_t<decltype(s)>,ApertureSource>) {
            auto launched=s.launch(scene,config.frequencyHz); states.insert(states.end(),launched.begin(),launched.end());
        } else states.push_back(s.launch(scene,config.frequencyHz));
    },source);
    return solve_states(scene,states,config);
}
Simulation solve(const Scene &scene, const Source &source, const SolverConfig &config) { return solve(scene,std::vector<Source>{source},config); }
} // namespace rt::beam
