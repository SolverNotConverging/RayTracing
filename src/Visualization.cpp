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

void visualize_rays(const std::vector<EMRay> &rays,
                    const std::vector<Rectangle> &rectangles,
                    const Vec3 &receiverPosition, double receiverRadius) {
    vtkNew<vtkPoints> pathPoints;
    vtkNew<vtkCellArray> pathLines;
    vtkNew<vtkPoints> polarizationPoints;
    vtkNew<vtkCellArray> polarizationLines;

    // Display scale only: the stored electric fields are unchanged.
    constexpr double polarizationScale = 0.12;
    constexpr int ellipseSegments = 64;

    for (const auto &ray : rays) {
        if (ray.path_.empty()) {
            continue;
        }

        // A polyline connects all saved samples, including intermediate steps.
        if (ray.path_.size() >= 2) {
            pathLines->InsertNextCell(static_cast<vtkIdType>(ray.path_.size()));
            for (const auto &sample : ray.path_) {
                pathLines->InsertCellPoint(pathPoints->InsertNextPoint(sample.position.data()));
            }
        }

        for (const auto &sample : ray.path_) {
            // Over one cycle, Re(E * exp(i*phase)) traces the polarization ellipse.
            // A linearly polarized field traces a line instead of an ellipse.
            const Vec3 realE = sample.E.real();
            const Vec3 imagE = sample.E.imag();
            const vtkIdType firstPoint = polarizationPoints->GetNumberOfPoints();
            polarizationLines->InsertNextCell(ellipseSegments + 1);
            for (int i = 0; i < ellipseSegments; ++i) {
                const double phase = 2.0 * PI * i / ellipseSegments;
                const Vec3 point = sample.position + polarizationScale *
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

    vtkNew<vtkPolyData> polarizationData;
    polarizationData->SetPoints(polarizationPoints);
    polarizationData->SetLines(polarizationLines);
    vtkNew<vtkPolyDataMapper> polarizationMapper;
    polarizationMapper->SetInputData(polarizationData);
    vtkNew<vtkActor> polarizationActor;
    polarizationActor->SetMapper(polarizationMapper);
    polarizationActor->GetProperty()->SetColor(1.0, 0.75, 0.2);
    polarizationActor->GetProperty()->SetLineWidth(2.0);

    // Use the same four corners and dimensions as the intersection geometry.
    vtkNew<vtkPoints> wallPoints;
    vtkNew<vtkCellArray> wallFaces;
    for (const Rectangle &wall : rectangles) {
        const Vec3 width = wall.halfWidth * wall.u;
        const Vec3 height = wall.halfHeight * wall.v;
        const Vec3 corners[] = {wall.center - width - height,
                                wall.center + width - height,
                                wall.center + width + height,
                                wall.center - width + height};
        wallFaces->InsertNextCell(4);
        for (const Vec3 &corner : corners) {
            wallFaces->InsertCellPoint(wallPoints->InsertNextPoint(corner.data()));
        }
    }
    vtkNew<vtkPolyData> wallData;
    wallData->SetPoints(wallPoints);
    wallData->SetPolys(wallFaces);
    vtkNew<vtkPolyDataMapper> wallMapper;
    wallMapper->SetInputData(wallData);
    vtkNew<vtkActor> wallActor;
    wallActor->SetMapper(wallMapper);
    wallActor->GetProperty()->SetColor(0.65, 0.7, 0.8);
    wallActor->GetProperty()->SetOpacity(0.25);
    wallActor->GetProperty()->EdgeVisibilityOn();
    wallActor->GetProperty()->SetEdgeColor(0.8, 0.85, 0.95);

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
    legend->SetInput("Cyan: ray paths | Gold: polarization (scaled)\n"
                     "Grey: PEC rectangles | Green: Rx reception sphere\n"
                     "Drag: rotate | Scroll: zoom | Middle drag: pan");
    legend->SetDisplayPosition(15, 15);
    legend->GetTextProperty()->SetFontSize(18);
    legend->GetTextProperty()->SetColor(1.0, 1.0, 1.0);

    vtkNew<vtkRenderer> renderer;
    renderer->SetBackground(0.08, 0.1, 0.14);
    renderer->AddActor(pathActor);
    renderer->AddActor(polarizationActor);
    renderer->AddActor(wallActor);
    renderer->AddActor(receiverActor);
    renderer->AddActor(axes);
    renderer->AddViewProp(legend);
    renderer->GetActiveCamera()->SetPosition(4.0, 3.0, 5.0);
    renderer->GetActiveCamera()->SetFocalPoint(0.0, 0.0, 0.0);
    renderer->GetActiveCamera()->SetViewUp(0.0, 0.0, 1.0);
    renderer->ResetCamera();

    vtkNew<vtkRenderWindow> window;
    window->SetWindowName("3D ray paths and polarization");
    window->SetSize(1100, 800);
    window->AddRenderer(renderer);

    vtkNew<vtkRenderWindowInteractor> interactor;
    vtkNew<vtkInteractorStyleTrackballCamera> style;
    interactor->SetRenderWindow(window);
    interactor->SetInteractorStyle(style);
    window->Render();
    interactor->Start();
}
