#include <raytracing/geometry.hpp>
#include <iostream>
#include <stdexcept>

using namespace rt;

namespace {
void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
void near(const Vec3 &a,const Vec3 &b,double tolerance,const char *message) {
    check((a-b).norm()<tolerance,message);
}
}

int main() try {
    const Vec3 X=Vec3::UnitX(),Y=Vec3::UnitY(),Z=Vec3::UnitZ();
    const EllipticCylinder ellipse{Vec3::Zero(),Z,X,1.0,0.6,2.0,false};
    validate_surface(ellipse);
    auto hit=intersect(Vec3::Zero(),X,ellipse);
    check(hit && std::abs(hit->distance-1)<1e-12,"Incorrect ellipse major semiaxis");
    near(hit->normal,X,1e-12,"Incorrect major-axis normal");
    hit=intersect(Vec3::Zero(),Y,ellipse);
    check(hit && std::abs(hit->distance-0.6)<1e-12,"Incorrect ellipse minor semiaxis");
    near(hit->normal,Y,1e-12,"Incorrect minor-axis normal");
    check(!intersect(Vec3::Zero(),Z,ellipse),"Open ellipse has a cap");
    auto closed=ellipse;closed.capped=true;
    check(intersect(Vec3::Zero(),Z,closed)->distance==2,"Ellipse cap location");
    check(!intersect(Vec3(0,0.7,3),-Z,closed),"Elliptical cap bounds");

    // Circular limit: independent established primitive must give the same hits.
    const Cylinder circle{Vec3::Zero(),Z,0.7,2,true};
    const EllipticCylinder circular{Vec3::Zero(),Z,X,0.7,0.7,2,true};
    for (int i=0;i<80;++i) {
        const Vec3 origin(0.1,0.15,0.3);
        const Vec3 direction=Vec3(std::cos(0.2*i),std::sin(0.2*i),0.1*i-4).normalized();
        const auto a=intersect(origin,direction,circle),b=intersect(origin,direction,circular);
        check(a && b && std::abs(a->distance-b->distance)<1e-11,"Circle-limit distance mismatch");
        near(a->normal,b->normal,1e-11,"Circle-limit normal mismatch");
    }

    // Curvature independently checked by differentiating the normalized normal.
    for (int i=0;i<20;++i) {
        const double t=0.1+0.29*i;
        const Vec3 p(std::cos(t),0.6*std::sin(t),0);
        const Vec3 tangent=Vec3(-std::sin(t),0.6*std::cos(t),0).normalized();
        auto normal=[](const Vec3 &q)->Vec3 {return Vec3(q.x(),q.y()/0.36,0).normalized();};
        const Vec3 n=normal(p);
        const auto h=curvature(ellipse,p,n);
        const Vec3 derivative=(normal(p+1e-6*tangent)-normal(p-1e-6*tangent))/2e-6;
        near(h*tangent,derivative,2e-8,"Ellipse curvature disagrees with normal derivative");
        near(h*Z,Vec3::Zero(),1e-12,"Cylinder has axial curvature");
        near(curvature(ellipse,p,-n)*tangent,-derivative,2e-8,"Inward curvature sign");
    }

    // The focus-to-focus optical property supplies a nontrivial path benchmark.
    // Each two-wall trajectory from F1 closes at F1 with transverse length 4a.
    const Vec3 focus(-0.8,0,0),otherFocus(0.8,0,0);
    for (double t : {0.17,0.49,1.1,1.9,2.7}) {
        Vec3 point=focus,direction(std::cos(t),std::sin(t),0);
        double length=0;
        for (int bounce=0;bounce<2;++bounce) {
            const auto event=intersect(point,direction,ellipse);
            check(event.has_value(),"Focus ray missed ellipse");
            point+=event->distance*direction;length+=event->distance;
            direction-=2*direction.dot(event->normal)*event->normal;
            if (bounce==0)
                near(direction,(otherFocus-point).normalized(),1e-10,"Ellipse focus reflection law");
        }
        near(direction,(focus-point).normalized(),1e-10,"Two-reflection focal closure");
        check(std::abs(length+(focus-point).norm()-4)<1e-10,"Focal round-trip length is not 4a");
    }

    const Annulus ring{Vec3(0,0,2),Z,0.4,1.0};
    validate_surface(ring);
    check(!intersect(Vec3::Zero(),Z,ring),"Annulus blocks its central opening");
    check(!intersect(Vec3(1.1,0,0),Z,ring),"Annulus outer bound");
    check(intersect(Vec3(0.7,0,0),Z,ring)->distance==2,"Annulus plane intersection");
    check(std::abs(edge_distance(ring,Vec3(0.7,0,2))-0.3)<1e-12,"Annulus edge distance");
    check(curvature(ring,Vec3(0.7,0,2),Z).norm()==0,"Annulus must be planar");
    bool rejected=false;
    try {validate_surface(Annulus{Vec3::Zero(),Z,1,0.5});} catch(const std::invalid_argument &) {rejected=true;}
    check(rejected,"Invalid annulus accepted");
    rejected=false;
    try {validate_surface(EllipticCylinder{Vec3::Zero(),Z,Z,1,0.6,2,false});}
    catch(const std::invalid_argument &) {rejected=true;}
    check(rejected,"Nonorthogonal ellipse axes accepted");
    std::cout << "PASS: elliptical intersections, circle limit, curvature, focal closure, and annular bounds\n";
    return 0;
} catch (const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
