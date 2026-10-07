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
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

void visualize_rays(const std::vector<EMRay> &rays) {
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

    vtkNew<vtkAxesActor> axes;
    axes->SetTotalLength(0.5, 0.5, 0.5);
    axes->SetXAxisLabelText("x [m]");
    axes->SetYAxisLabelText("y [m]");
    axes->SetZAxisLabelText("z [m]");

    vtkNew<vtkTextActor> legend;
    legend->SetInput("Cyan: ray paths\nGold: polarization at each sample (scaled for display)\n"
                     "Drag: rotate | Scroll: zoom | Middle drag: pan");
    legend->SetDisplayPosition(15, 15);
    legend->GetTextProperty()->SetFontSize(18);
    legend->GetTextProperty()->SetColor(1.0, 1.0, 1.0);

    vtkNew<vtkRenderer> renderer;
    renderer->SetBackground(0.08, 0.1, 0.14);
    renderer->AddActor(pathActor);
    renderer->AddActor(polarizationActor);
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
