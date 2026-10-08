#include "Visualization.hpp"

#include <cmath>
#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSphereSource.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <type_traits>
#include <vtkSmartPointer.h>

namespace {
    // These meshes are only for display; ray intersections use Surface.cpp.
    vtkSmartPointer<vtkPolyData> surface_mesh(const Surface &surface) {
        return std::visit([](const auto &shape) -> vtkSmartPointer<vtkPolyData> {
            using Shape = std::decay_t<decltype(shape)>;
            if constexpr (std::is_same_v<Shape, Sphere>) {
                vtkNew<vtkSphereSource> sphere;
                sphere->SetCenter(shape.center_.data());
                sphere->SetRadius(shape.radius_);
                sphere->SetThetaResolution(64);
                sphere->SetPhiResolution(32);
                sphere->Update();
                return sphere->GetOutput();
            } else {
                vtkNew<vtkPoints> points;
                vtkNew<vtkCellArray> faces;
                auto polygon = [&](const std::vector<Vec3> &vertices) {
                    faces->InsertNextCell(static_cast<vtkIdType>(vertices.size()));
                    for (const Vec3 &vertex: vertices)
                        faces->InsertCellPoint(points->InsertNextPoint(vertex.data()));
                };
                if constexpr (std::is_same_v<Shape, Rectangle>) {
                    const Vec3 u = shape.halfWidth_ * shape.u_, v = shape.halfHeight_ * shape.v_;
                    polygon({
                        shape.center_ - u - v, shape.center_ + u - v,
                        shape.center_ + u + v, shape.center_ - u + v
                    });
                } else if constexpr (std::is_same_v<Shape, Triangle>) {
                    polygon({shape.a_, shape.b_, shape.c_});
                } else {
                    constexpr int segments = 64;
                    const Vec3 axis = [&]() -> Vec3 {
                        if constexpr (std::is_same_v<Shape, Disk>) return shape.normal_;
                        else return shape.axis_;
                    }();
                    const Vec3 u = axis.unitOrthogonal();
                    const Vec3 v = axis.cross(u);
                    std::vector<Vec3> ring;
                    for (int i = 0; i < segments; ++i) {
                        const double angle = 2.0 * PI * i / segments;
                        ring.push_back(shape.center_ + shape.radius_ *
                                       (std::cos(angle) * u + std::sin(angle) * v));
                    }
                    if constexpr (std::is_same_v<Shape, Disk>) {
                        polygon(ring);
                    } else {
                        const Vec3 offset = shape.halfLength_ * axis;
                        for (int i = 0; i < segments; ++i) {
                            const int j = (i + 1) % segments;
                            polygon({
                                ring[i] - offset, ring[j] - offset,
                                ring[j] + offset, ring[i] + offset
                            });
                        }
                        if (shape.capped_) {
                            std::vector<Vec3> bottom, top;
                            for (int i = 0; i < segments; ++i) {
                                bottom.push_back(ring[segments - 1 - i] - offset);
                                top.push_back(ring[i] + offset);
                            }
                            polygon(bottom);
                            polygon(top);
                        }
                    }
                }
                auto data = vtkSmartPointer<vtkPolyData>::New();
                data->SetPoints(points);
                data->SetPolys(faces);
                return data;
            }
        }, surface);
    }
}

static void show_scene(const std::vector<std::vector<Vec3>> &paths,
                    const std::vector<EMRay> &rays,
                    const std::vector<Surface> &surfaces,
                    const Vec3 &receiverPosition, double receiverRadius,
                    const std::vector<std::vector<Vec3>> &unresolvedPaths = {}) {
    vtkNew<vtkPoints> pathPoints;
    vtkNew<vtkCellArray> pathLines;
    vtkNew<vtkPoints> polarizationPoints;
    vtkNew<vtkCellArray> polarizationLines;

    // Display scale only: the stored electric fields are unchanged.
    constexpr double polarizationScale = 0.12;
    constexpr int ellipseSegments = 64;

    for (const auto &path : paths) {
        if (path.size() < 2) continue;
        pathLines->InsertNextCell(static_cast<vtkIdType>(path.size()));
        for (const Vec3 &point : path)
            pathLines->InsertCellPoint(pathPoints->InsertNextPoint(point.data()));
    }

    for (const auto &ray: rays) {
        if (ray.path_.empty()) {
            continue;
        }

        for (const auto &sample: ray.path_) {
            // Over one cycle, Re(E * exp(i*phase)) traces the polarization ellipse.
            // A linearly polarized field traces a line instead of an ellipse.
            const Vec3 realE = sample.E_.real();
            const Vec3 imagE = sample.E_.imag();
            const vtkIdType firstPoint = polarizationPoints->GetNumberOfPoints();
            polarizationLines->InsertNextCell(ellipseSegments + 1);
            for (int i = 0; i < ellipseSegments; ++i) {
                const double phase = 2.0 * PI * i / ellipseSegments;
                const Vec3 point = sample.position_ + polarizationScale *
                                   (realE * std::cos(phase) - imagE * std::sin(phase));
                polarizationLines->InsertCellPoint(polarizationPoints->InsertNextPoint(point.data()));
            }
            polarizationLines->InsertCellPoint(firstPoint);
        }
    }

    vtkNew<vtkPolyData> pathData;
    pathData->SetPoints(pathPoints);
    pathData->SetLines(pathLines);
    vtkNew<vtkPolyDataMapper> pathMapper;
    pathMapper->SetInputData(pathData);
    vtkNew<vtkActor> pathActor;
    pathActor->SetMapper(pathMapper);
    pathActor->GetProperty()->SetColor(0.2, 0.8, 1.0);
    pathActor->GetProperty()->SetLineWidth(2.0);

    vtkNew<vtkPoints> unresolvedPoints;
    vtkNew<vtkCellArray> unresolvedLines;
    for (const auto &path : unresolvedPaths) {
        if (path.size() < 2) continue;
        unresolvedLines->InsertNextCell(static_cast<vtkIdType>(path.size()));
        for (const Vec3 &point : path)
            unresolvedLines->InsertCellPoint(unresolvedPoints->InsertNextPoint(point.data()));
    }
    vtkNew<vtkPolyData> unresolvedData;
    unresolvedData->SetPoints(unresolvedPoints);
    unresolvedData->SetLines(unresolvedLines);
    vtkNew<vtkPolyDataMapper> unresolvedMapper;
    unresolvedMapper->SetInputData(unresolvedData);
    vtkNew<vtkActor> unresolvedActor;
    unresolvedActor->SetMapper(unresolvedMapper);
    unresolvedActor->GetProperty()->SetColor(1.0, 0.45, 0.1);
    unresolvedActor->GetProperty()->SetLineWidth(2.0);

    vtkNew<vtkPolyData> polarizationData;
    polarizationData->SetPoints(polarizationPoints);
    polarizationData->SetLines(polarizationLines);
    vtkNew<vtkPolyDataMapper> polarizationMapper;
    polarizationMapper->SetInputData(polarizationData);
    vtkNew<vtkActor> polarizationActor;
    polarizationActor->SetMapper(polarizationMapper);
    polarizationActor->GetProperty()->SetColor(1.0, 0.75, 0.2);
    polarizationActor->GetProperty()->SetLineWidth(2.0);

    vtkNew<vtkSphereSource> receiverSphere;
    receiverSphere->SetCenter(receiverPosition.data());
    receiverSphere->SetRadius(receiverRadius);
    receiverSphere->SetThetaResolution(32);
    receiverSphere->SetPhiResolution(24);
    vtkNew<vtkPolyDataMapper> receiverMapper;
    receiverMapper->SetInputConnection(receiverSphere->GetOutputPort());
    vtkNew<vtkActor> receiverActor;
    receiverActor->SetMapper(receiverMapper);
    receiverActor->GetProperty()->SetColor(0.2, 1.0, 0.35);
    receiverActor->GetProperty()->SetOpacity(0.5);

    vtkNew<vtkAxesActor> axes;
    axes->SetTotalLength(0.5, 0.5, 0.5);
    axes->SetXAxisLabelText("x [m]");
    axes->SetYAxisLabelText("y [m]");
    axes->SetZAxisLabelText("z [m]");

    vtkNew<vtkTextActor> legend;
    legend->SetInput(rays.empty() ?
        "Cyan: refined paths (geometry only)\n"
        "Orange: unresolved corner candidates (coarse paths)\n"
        "Grey: PEC surfaces | Green: Rx reception sphere\n"
        "Drag: rotate | Scroll: zoom | Middle drag: pan" :
        "Cyan: ray paths | Gold: polarization (scaled)\n"
        "Grey: PEC surfaces | Green: Rx reception sphere\n"
        "Drag: rotate | Scroll: zoom | Middle drag: pan");
    legend->SetDisplayPosition(15, 15);
    legend->GetTextProperty()->SetFontSize(18);
    legend->GetTextProperty()->SetColor(1.0, 1.0, 1.0);

    vtkNew<vtkRenderer> renderer;
    renderer->SetBackground(0.08, 0.1, 0.14);
    renderer->AddActor(pathActor);
    renderer->AddActor(unresolvedActor);
    renderer->AddActor(polarizationActor);
    for (const Surface &surface: surfaces) {
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(surface_mesh(surface));
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(0.65, 0.7, 0.8);
        actor->GetProperty()->SetOpacity(0.25);
        if (std::holds_alternative<Rectangle>(surface) || std::holds_alternative<Triangle>(surface)) {
            actor->GetProperty()->EdgeVisibilityOn();
            actor->GetProperty()->SetEdgeColor(0.8, 0.85, 0.95);
        }
        renderer->AddActor(actor);
    }
    renderer->AddActor(receiverActor);
    renderer->AddActor(axes);
    renderer->AddViewProp(legend);
    renderer->GetActiveCamera()->SetPosition(4.0, 3.0, 5.0);
    renderer->GetActiveCamera()->SetFocalPoint(0.0, 0.0, 0.0);
    renderer->GetActiveCamera()->SetViewUp(0.0, 0.0, 1.0);
    renderer->ResetCamera();

    vtkNew<vtkRenderWindow> window;
    window->SetWindowName(rays.empty() ? "Refined ray geometry" : "3D ray paths and polarization");
    window->SetSize(1100, 800);
    window->AddRenderer(renderer);

    vtkNew<vtkRenderWindowInteractor> interactor;
    vtkNew<vtkInteractorStyleTrackballCamera> style;
    interactor->SetRenderWindow(window);
    interactor->SetInteractorStyle(style);
    window->Render();
    interactor->Start();
}

void visualize_rays(const std::vector<EMRay> &rays, const std::vector<Surface> &surfaces,
                    const Vec3 &receiverPosition, double receiverRadius) {
    std::vector<std::vector<Vec3>> paths;
    for (const auto &ray : rays) {
        std::vector<Vec3> points;
        for (const auto &sample : ray.path_) points.push_back(sample.position_);
        paths.push_back(std::move(points));
    }
    show_scene(paths, rays, surfaces, receiverPosition, receiverRadius);
}

void visualize_paths(const std::vector<std::vector<Vec3>> &paths, const std::vector<Surface> &surfaces,
                     const Vec3 &receiverPosition, double receiverRadius,
                     const std::vector<std::vector<Vec3>> &unresolvedPaths) {
    show_scene(paths, {}, surfaces, receiverPosition, receiverRadius, unresolvedPaths);
}
