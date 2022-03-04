#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkContourFilter.h>
#include <vtkSmartPointer.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <stdio.h>
#include <iostream>
#include <vtkCamera.h>

int main(int argc, char* argv[])
{
    vtkDataSetReader* reader = vtkDataSetReader::New();
    reader->SetFileName("noise.vtk");
    reader->Update();

    vtkDataSetMapper* mapper = vtkDataSetMapper::New();
    mapper->SetInputData(reader->GetOutput());

    vtkLookupTable* lut = vtkLookupTable::New();

    //Part D
    int n = 256;
    lut->SetNumberOfTableValues(n);
    for (int i = 0; i < n; i++) {
        lut->SetTableValue(n - 1 - i, (i / static_cast<float>(n)), 0.0, 1 - (i / static_cast<float>(n)), 1.0);
    }
    lut->Build();

    mapper->SetLookupTable(lut);
    mapper->SetScalarRange(1, 6);

    vtkActor* actor = vtkActor::New();
    actor->SetMapper(mapper);
    //end Part D


    //Part B
    vtkSmartPointer<vtkContourFilter> contour_filter = vtkSmartPointer<vtkContourFilter>::New();
    contour_filter->SetNumberOfContours(2);
    contour_filter->SetInputData(reader->GetOutput());
    contour_filter->SetValue(1, 1);

    vtkSmartPointer<vtkPolyDataMapper> contour_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    contour_mapper->SetInputConnection(contour_filter->GetOutputPort());
    contour_mapper->SetLookupTable(lut);
    contour_mapper->SetScalarRange(1, 6);

    vtkSmartPointer<vtkActor> contour_actor = vtkSmartPointer<vtkActor>::New();
    contour_actor->SetMapper(contour_mapper);
    //end Part B

    //Part C
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(0.0, 0.0, 0.0);
    plane->SetNormal(0.0, 0.0, 1.0);


    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputConnection(reader->GetOutputPort());
    cutter->SetCutFunction(plane);
    cutter->Update();

    vtkSmartPointer<vtkPolyDataMapper> cut_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cut_mapper->SetInputConnection(cutter->GetOutputPort());
    cut_mapper->SetLookupTable(lut);
    cut_mapper->SetScalarRange(1, 6);

    vtkSmartPointer<vtkActor> cut_actor = vtkSmartPointer<vtkActor>::New();
    cut_actor->SetMapper(cut_mapper);
    //end Part C


    vtkRenderer* ren = vtkRenderer::New();
    //ren->AddActor(actor);
    ren->AddActor(cut_actor); //Part C/D/E
    ren->SetViewport(0.0, 0.0, 0.5, 1.0);

    //Part F
    vtkRenderer* ren2 = vtkRenderer::New();
    ren2->AddActor(contour_actor); //Part B/E
    ren2->SetViewport(0.5, 0.0, 1.0, 1.0);
    //end Part F
    vtkRenderWindow* renwin = vtkRenderWindow::New();
    renwin->AddRenderer(ren);
    renwin->AddRenderer(ren2);//Part F
    renwin->SetSize(1536, 768); //Part A

    //vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    //iren->SetRenderWindow(renwin);
    renwin->Render();
    ren2->ResetCamera();
    vtkSmartPointer<vtkCamera> c2 = ren2->GetActiveCamera();
    c2->SetPosition(0.0,0.0,67.0);

    int i = 0;
    for (i; i < 500; i++) {
        contour_filter->SetValue(1, 1 + (i * 0.01));
        contour_filter->Update();        
        renwin->Render();
    }
    
    //iren->Start();
}
