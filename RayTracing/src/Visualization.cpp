#include "Visualization.hpp"
#include "ResultIO.hpp"

#include <cmath>
#include <stdexcept>
#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>
#include <vtkCellArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkWindowToImageFilter.h>

#include <type_traits>
#include <vtkSmartPointer.h>

namespace {
    // These meshes are only for display; ray intersections use Surface.cpp.
    vtkSmartPointer<vtkPolyData> surface_mesh(const Surface &surface) {
        return std::visit(
            [](const auto &shape) -> vtkSmartPointer<vtkPolyData> {
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
                        const Vec3 u = shape.halfWidth_ * shape.u_,
                                v = shape.halfHeight_ * shape.v_;
                        polygon({
                            shape.center_ - u - v, shape.center_ + u - v,
                            shape.center_ + u + v, shape.center_ - u + v
                        });
                    } else if constexpr (std::is_same_v<Shape, Triangle>) {
                        polygon({shape.a_, shape.b_, shape.c_});
                    } else {
                        constexpr int segments = 64;
                        const Vec3 axis = [&]() -> Vec3 {
                            if constexpr (std::is_same_v<Shape, Disk>)
                                return shape.normal_;
                            else
                                return shape.axis_;
                        }();
                        const Vec3 u = axis.unitOrthogonal();
                        const Vec3 v = axis.cross(u);
                        std::vector<Vec3> ring;
                        for (int i = 0; i < segments; ++i) {
                            const double angle = 2.0 * PI * i / segments;
                            ring.push_back(shape.center_ +
                                           shape.radius_ *
                                           (std::cos(angle) * u + std::sin(angle) * v));
                        }
                        if constexpr (std::is_same_v<Shape, Disk>) {
                            polygon(ring);
                        } else {
                            const Vec3 offset = shape.halfLength_ * axis;
                            for (int i = 0; i < segments; ++i) {
                                const int j = (i + 1) % segments;
                                polygon({
                                    ring[i] - offset, ring[j] - offset, ring[j] + offset,
                                    ring[i] + offset
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
            },
            surface);
    }
} // namespace

static vtkSmartPointer<vtkPolyData>
antenna_mesh(const rt::Antenna &antenna, const Vec3 &position, double frequency,
             const rt::Medium &medium, double scale) {
    constexpr int thetaCount = 61, phiCount = 120;
    std::vector<Vec3> directions;
    std::vector<double> amplitudes;
    double peak = 0;
    double thetaStart = 0, thetaEnd = 180, phiStart = 0, phiSpan = 360;
    bool periodic = true;
    if (antenna.pattern) {
        thetaStart = antenna.pattern->thetaDegrees.front();
        thetaEnd = antenna.pattern->thetaDegrees.back();
        const auto &axis = antenna.pattern->phiDegrees;
        phiStart = axis.front();
        phiSpan = axis.back() - axis.front();
        periodic = 360 - phiSpan <= 1.01 * (axis.back() - axis[axis.size() - 2]);
        if (periodic)
            phiSpan = 360;
    }
    for (int t = 0; t < thetaCount; ++t)
        for (int p = 0; p < phiCount; ++p) {
            const double theta =
                    (thetaStart + (thetaEnd - thetaStart) * t / (thetaCount - 1)) * PI /
                    180;
            const double phi =
                    (phiStart + phiSpan * p / (periodic ? phiCount : phiCount - 1)) * PI /
                    180;
            const Vec3 d = antenna.orientation * Vec3(std::sin(theta) * std::cos(phi),
                                                      std::sin(theta) * std::sin(phi),
                                                      std::cos(theta));
            const double amplitude = antenna.farfield(d, frequency, medium).norm();
            directions.push_back(d);
            amplitudes.push_back(amplitude);
            peak = std::max(peak, amplitude);
        }
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> faces;
    for (std::size_t i = 0; i < directions.size(); ++i) {
        const Vec3 p = position + (peak > 0 ? scale * amplitudes[i] / peak : 0) *
                       directions[i];
        points->InsertNextPoint(p.data());
    }
    for (int t = 0; t < thetaCount - 1; ++t)
        for (int p = 0; p < (periodic ? phiCount : phiCount - 1); ++p) {
            const int q = (p + 1) % phiCount;
            const vtkIdType vertices[4] = {
                t * phiCount + p, t * phiCount + q,
                (t + 1) * phiCount + q,
                (t + 1) * phiCount + p
            };
            faces->InsertNextCell(4, vertices);
        }
    auto mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->SetPoints(points);
    mesh->SetPolys(faces);
    return mesh;
}

static void
show_scene(const std::vector<std::vector<Vec3> > &paths,
           const std::vector<EMRay> &rays, const std::vector<Surface> &surfaces,
           const Vec3 &receiverPosition, double receiverRadius,
           const std::vector<std::vector<Vec3> > &unresolvedPaths = {},
           const rt::SimulationResult *simulation = nullptr,
           const rt::ViewOptions &view = {}) {
    vtkNew<vtkPoints> pathPoints;
    vtkNew<vtkCellArray> pathLines;
    vtkNew<vtkPoints> polarizationPoints;
    vtkNew<vtkCellArray> polarizationLines;

    // Display scale only: the stored electric fields are unchanged.
    constexpr double polarizationScale = 0.12;
    constexpr int ellipseSegments = 64;

    for (const auto &path: paths) {
        if (path.size() < 2)
            continue;
        pathLines->InsertNextCell(static_cast<vtkIdType>(path.size()));
        for (const Vec3 &point: path)
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
                const Vec3 point =
                        sample.position_ + polarizationScale * (realE * std::cos(phase) -
                                                                imagE * std::sin(phase));
                polarizationLines->InsertCellPoint(
                    polarizationPoints->InsertNextPoint(point.data()));
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
    for (const auto &path: unresolvedPaths) {
        if (path.size() < 2)
            continue;
        unresolvedLines->InsertNextCell(static_cast<vtkIdType>(path.size()));
        for (const Vec3 &point: path)
            unresolvedLines->InsertCellPoint(
                unresolvedPoints->InsertNextPoint(point.data()));
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
    if (simulation &&
        simulation->receiver.antenna.kind != rt::AntennaKind::Isotropic)
        receiverMapper->SetInputData(
            antenna_mesh(simulation->receiver.antenna, receiverPosition,
                         simulation->settings.frequencyHz,
                         simulation->settings.medium, view.patternScale));
    else
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
    for (auto caption:
         {
             axes->GetXAxisCaptionActor2D(), axes->GetYAxisCaptionActor2D(),
             axes->GetZAxisCaptionActor2D()
         }) {
        caption->GetTextActor()->SetTextScaleModeToNone();
        caption->GetCaptionTextProperty()->SetFontSize(14);
    }

    vtkNew<vtkTextActor> legend;
    legend->SetInput(rays.empty()
                         ? "Cyan: refined paths (geometry only)\n"
                         "Orange: unresolved corner candidates (coarse paths)\n"
                         "Grey: PEC surfaces | Green: Rx reception sphere\n"
                         "Drag: rotate | Scroll: zoom | Middle drag: pan"
                         : "Cyan: ray paths | Gold: polarization (scaled)\n"
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
        if (std::holds_alternative<Rectangle>(surface) ||
            std::holds_alternative<Triangle>(surface)) {
            actor->GetProperty()->EdgeVisibilityOn();
            actor->GetProperty()->SetEdgeColor(0.8, 0.85, 0.95);
        }
        renderer->AddActor(actor);
    }
    renderer->AddActor(receiverActor);
    if (simulation) {
        vtkNew<vtkPolyDataMapper> txMapper;
        vtkNew<vtkSphereSource> txSphere;
        if (simulation->transmitter.antenna.kind != rt::AntennaKind::Isotropic)
            txMapper->SetInputData(antenna_mesh(
                simulation->transmitter.antenna, simulation->transmitter.position,
                simulation->settings.frequencyHz, simulation->settings.medium,
                view.patternScale));
        else {
            txSphere->SetCenter(simulation->transmitter.position.data());
            txSphere->SetRadius(receiverRadius);
            txSphere->SetThetaResolution(32);
            txSphere->SetPhiResolution(24);
            txMapper->SetInputConnection(txSphere->GetOutputPort());
        }
        vtkNew<vtkActor> txActor;
        txActor->SetMapper(txMapper);
        txActor->GetProperty()->SetColor(1, 0.3, 0.25);
        txActor->GetProperty()->SetOpacity(0.75);
        renderer->AddActor(txActor);
        legend->SetInput("Cyan: refined paths | Gold: incident field polarization\n"
            "Red: Tx | Green: Rx (normalized antenna-pattern shape)\n"
            "Orange: unresolved candidates | Grey: PEC geometry\n"
            "Drag: rotate | Scroll: zoom | Middle drag: pan");
    }
    renderer->AddActor(axes);
    renderer->AddViewProp(legend);
    renderer->GetActiveCamera()->SetPosition(4.0, 3.0, 5.0);
    renderer->GetActiveCamera()->SetFocalPoint(0.0, 0.0, 0.0);
    renderer->GetActiveCamera()->SetViewUp(0.0, 0.0, 1.0);
    renderer->ResetCamera();

    vtkNew<vtkRenderWindow> window;
    window->SetWindowName(rays.empty()
                              ? "Refined ray geometry"
                              : "3D ray paths and polarization");
    window->SetSize(1100, 800);
    window->AddRenderer(renderer);
    if (!view.interactive)
        window->SetOffScreenRendering(1);

    vtkNew<vtkRenderWindowInteractor> interactor;
    vtkNew<vtkInteractorStyleTrackballCamera> style;
    interactor->SetRenderWindow(window);
    interactor->SetInteractorStyle(style);
    window->Render();
    if (!view.screenshot.empty()) {
        if (!view.screenshot.parent_path().empty())
            std::filesystem::create_directories(view.screenshot.parent_path());
        vtkNew<vtkWindowToImageFilter> capture;
        capture->SetInput(window);
        capture->ReadFrontBufferOff();
        capture->Update();
        vtkNew<vtkPNGWriter> writer;
        writer->SetFileName(view.screenshot.string().c_str());
        writer->SetInputConnection(capture->GetOutputPort());
        writer->Write();
        if (!std::filesystem::exists(view.screenshot))
            throw std::runtime_error("Failed exporting VTK screenshot");
    }
    if (view.interactive)
        interactor->Start();
}

namespace rt {
    void visualize(const SimulationResult &result, const ViewOptions &options) {
        if (!std::isfinite(options.patternScale) || options.patternScale <= 0)
            throw std::invalid_argument("Pattern display scale must be positive");
        std::vector<std::vector<Vec3> > paths, unresolved;
        std::vector<EMRay> fields;
        auto indices = options.rayIndices;
        if (indices.empty())
            for (std::size_t i = 0; i < result.rays.size(); ++i)
                indices.push_back(i);
        for (auto index: indices) {
            const auto &ray = result.rays.at(index);
            const auto &geometry = ray.refinement.geometry_;
            std::vector<Vec3> vertices{result.transmitter.position};
            for (const auto &hit: geometry.reflections_)
                vertices.push_back(hit.position_);
            if (geometry.receiver_)
                vertices.push_back(geometry.receiver_->closestPoint_);
            paths.push_back(std::move(vertices));
            if (ray.field) {
                EMRay display(result.transmitter.position,
                              ray.refinement.launchDirection_,
                              launch_polarization(ray.refinement.launchDirection_),
                              C / result.settings.frequencyHz);
                display.path_ = {
                    {
                        result.receiver.position, ray.field->receiverField_,
                        ray.field->opticalPath_
                    }
                };
                fields.push_back(std::move(display));
            }
        }
        for (const auto &candidate: result.candidates)
            if (candidate.refinement.status_ == RefinementStatus::UnresolvedCorner)
                unresolved.push_back(candidate.coarseVertices);
        show_scene(paths, fields, result.scene.surfaces(), result.receiver.position,
                   result.settings.receptionRadius, unresolved, &result, options);
    }

    void visualize(const SimulationResult &result, std::size_t index,
                   const ViewOptions &options) {
        auto selected = options;
        selected.rayIndices = {index};
        visualize(result, selected);
    }

    void visualize_h5(const std::filesystem::path &file,
                      const ViewOptions &options) {
        visualize(load_h5(file), options);
    }
} // namespace rt

void visualize_rays(const std::vector<EMRay> &rays,
                    const std::vector<Surface> &surfaces,
                    const Vec3 &receiverPosition, double receiverRadius) {
    std::vector<std::vector<Vec3> > paths;
    for (const auto &ray: rays) {
        std::vector<Vec3> points;
        for (const auto &sample: ray.path_)
            points.push_back(sample.position_);
        paths.push_back(std::move(points));
    }
    show_scene(paths, rays, surfaces, receiverPosition, receiverRadius);
}

void visualize_paths(const std::vector<std::vector<Vec3> > &paths,
                     const std::vector<Surface> &surfaces,
                     const Vec3 &receiverPosition, double receiverRadius,
                     const std::vector<std::vector<Vec3> > &unresolvedPaths) {
    show_scene(paths, {}, surfaces, receiverPosition, receiverRadius,
               unresolvedPaths);
}
