#include "visualization.hpp"
#include <RayTracing/Constants.hpp>
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSphereSource.h>
#include <vtkWindowToImageFilter.h>
#include <cmath>
#include <stdexcept>

namespace {
struct Mesh {
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> faces;
    vtkIdType point(const Vec3 &p) { return points->InsertNextPoint(p.x(), p.y(), p.z()); }
    void triangle(const Vec3 &a, const Vec3 &b, const Vec3 &c) {
        vtkIdType ids[] = {point(a), point(b), point(c)};
        faces->InsertNextCell(3, ids);
    }
    void quad(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d) {
        triangle(a, b, c); triangle(a, c, d);
    }
};
constexpr int resolution = 40;
constexpr double tau = 2 * rt::constants::pi;
void surface_mesh(Mesh &mesh, const Rectangle &s) {
    const Vec3 u = s.halfWidth_ * s.u_, v = s.halfHeight_ * s.v_;
    mesh.quad(s.center_ - u - v, s.center_ + u - v, s.center_ + u + v, s.center_ - u + v);
}
void surface_mesh(Mesh &mesh, const Triangle &s) { mesh.triangle(s.a_, s.b_, s.c_); }
void surface_mesh(Mesh &mesh, const Disk &s) {
    const Vec3 u = s.normal_.unitOrthogonal(), v = s.normal_.cross(u);
    for (int i = 0; i < resolution; ++i) {
        const double a = tau * i / resolution, b = tau * (i + 1) / resolution;
        mesh.triangle(s.center_, s.center_ + s.radius_ * (std::cos(a) * u + std::sin(a) * v),
                      s.center_ + s.radius_ * (std::cos(b) * u + std::sin(b) * v));
    }
}
void surface_mesh(Mesh &mesh, const Cylinder &s) {
    const Vec3 u = s.axis_.unitOrthogonal(), v = s.axis_.cross(u), dz = s.halfLength_ * s.axis_;
    for (int i = 0; i < resolution; ++i) {
        const double a = tau * i / resolution, b = tau * (i + 1) / resolution;
        const Vec3 p = s.center_ + s.radius_ * (std::cos(a) * u + std::sin(a) * v);
        const Vec3 q = s.center_ + s.radius_ * (std::cos(b) * u + std::sin(b) * v);
        mesh.quad(p - dz, q - dz, q + dz, p + dz);
        if (s.capped_) {
            mesh.triangle(s.center_ - dz, q - dz, p - dz);
            mesh.triangle(s.center_ + dz, p + dz, q + dz);
        }
    }
}
void surface_mesh(Mesh &mesh, const Sphere &s) {
    auto point = [&](int latitude, int longitude) -> Vec3 {
        const double theta = rt::constants::pi * latitude / (resolution / 2);
        const double phi = tau * longitude / resolution;
        return s.center_ + s.radius_ * Vec3(std::sin(theta) * std::cos(phi),
                                          std::sin(theta) * std::sin(phi), std::cos(theta));
    };
    for (int j = 0; j < resolution / 2; ++j)
        for (int i = 0; i < resolution; ++i)
            mesh.quad(point(j, i), point(j, i + 1), point(j + 1, i + 1), point(j + 1, i));
}
void add_actor(vtkRenderer *renderer, vtkPolyData *data, double r, double g, double b,
               double opacity = 1) {
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(data);
    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(r, g, b);
    actor->GetProperty()->SetOpacity(opacity);
    actor->GetProperty()->SetLineWidth(2);
    renderer->AddActor(actor);
}
void marker(vtkRenderer *renderer, const Vec3 &p, double r, double g, double b) {
    vtkNew<vtkSphereSource> sphere;
    sphere->SetCenter(p.x(), p.y(), p.z());
    sphere->SetRadius(0.07);
    sphere->Update();
    add_actor(renderer, sphere->GetOutput(), r, g, b);
}
}

void view_rays(const rt::SimulationResult &result, bool interactive,
               const std::filesystem::path &screenshot) {
    vtkNew<vtkRenderer> renderer;
    renderer->SetBackground(0.96, 0.97, 1);
    Mesh scene;
    for (const auto &surface : result.scene.surfaces())
        std::visit([&](const auto &s) { surface_mesh(scene, s); }, surface);
    vtkNew<vtkPolyData> geometry;
    geometry->SetPoints(scene.points);
    geometry->SetPolys(scene.faces);
    add_actor(renderer, geometry, 0.55, 0.60, 0.67, 0.22);

    // Cyan: valid field paths. Orange: refined geometry with invalid spreading.
    // Each polyline ends at the refined receiver point, not the capture sphere.
    for (bool valid : {false, true}) {
        vtkNew<vtkPoints> points;
        vtkNew<vtkCellArray> lines;
        for (const auto &ray : result.rays) {
            if (ray.field.has_value() != valid) continue;
            const auto &path = ray.refinement.geometry_;
            if (!path.receiver_) continue;
            lines->InsertNextCell(static_cast<vtkIdType>(path.reflections_.size() + 2));
            auto append = [&](const Vec3 &p) {
                lines->InsertCellPoint(points->InsertNextPoint(p.x(), p.y(), p.z()));
            };
            append(result.transmitter.position);
            for (const auto &hit : path.reflections_) append(hit.position_);
            append(path.receiver_->closestPoint_);
        }
        vtkNew<vtkPolyData> paths;
        paths->SetPoints(points);
        paths->SetLines(lines);
        if (valid) add_actor(renderer, paths, 0.02, 0.50, 0.72);
        else add_actor(renderer, paths, 0.95, 0.50, 0.05);
    }
    marker(renderer, result.transmitter.position, 0.85, 0.1, 0.1);
    marker(renderer, result.receiver.position, 0.1, 0.7, 0.2);
    vtkNew<vtkRenderWindow> window;
    window->AddRenderer(renderer);
    window->SetSize(1200, 850);
    window->SetWindowName("Point rays: Tx red, Rx green, valid cyan, invalid field orange");
    window->SetOffScreenRendering(!interactive);
    renderer->ResetCamera();
    window->Render();
    if (!screenshot.empty()) {
        vtkNew<vtkWindowToImageFilter> capture;
        capture->SetInput(window);
        capture->ReadFrontBufferOff();
        capture->Update();
        vtkNew<vtkPNGWriter> writer;
        writer->SetFileName(screenshot.string().c_str());
        writer->SetInputConnection(capture->GetOutputPort());
        writer->Write();
        if (writer->GetErrorCode() != 0) throw std::runtime_error("VTK could not write screenshot");
    }
    if (interactive) {
        vtkNew<vtkRenderWindowInteractor> interactor;
        interactor->SetRenderWindow(window);
        interactor->Initialize();
        interactor->Start();
    }
}
