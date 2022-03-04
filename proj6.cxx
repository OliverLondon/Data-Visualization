/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"					//
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"				//
#include "vtkRenderWindow.h"			//
#include "vtkRenderWindowInteractor.h"	//
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>			//
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>			//
#include <vtkRenderer.h>				//
#include <vtkActor.h>					//
#include <vtkRenderWindow.h> 			//
#include <vtkRenderWindowInteractor.h>	//
#include <vtkSmartPointer.h>

#include <cmath>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    // idx[0] = pointId%dims[0];
    // idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}


class TriangleList
{
   public:
                   TriangleList() { maxTriangles = 1000000; triangleIdx = 0; pts = new float[9*maxTriangles]; };
     virtual      ~TriangleList() { delete [] pts; };

     void          AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxTriangles;
     int           triangleIdx;
};

void
TriangleList::AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3)
{
    pts[9*triangleIdx+0] = X1;
    pts[9*triangleIdx+1] = Y1;
    pts[9*triangleIdx+2] = Z1;
    pts[9*triangleIdx+3] = X2;
    pts[9*triangleIdx+4] = Y2;
    pts[9*triangleIdx+5] = Z2;
    pts[9*triangleIdx+6] = X3;
    pts[9*triangleIdx+7] = Y3;
    pts[9*triangleIdx+8] = Z3;
    triangleIdx++;
}

vtkPolyData *
TriangleList::MakePolyData(void)
{
    int ntriangles = triangleIdx;
    int numPoints = 3*(ntriangles);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *tris = vtkCellArray::New();
    tris->EstimateSize(numPoints,4);
    for (int i = 0 ; i < ntriangles ; i++)
    {
        double pt[3];
        pt[0] = pts[9*i];
        pt[1] = pts[9*i+1];
        pt[2] = pts[9*i+2];
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[9*i+3];
        pt[1] = pts[9*i+4];
        pt[2] = pts[9*i+5];
        vtk_pts->SetPoint(ptIdx+1, pt);
        pt[0] = pts[9*i+6];
        pt[1] = pts[9*i+7];
        pt[2] = pts[9*i+8];
        vtk_pts->SetPoint(ptIdx+2, pt);
        vtkIdType ids[3] = { ptIdx, ptIdx+1, ptIdx+2 };
        tris->InsertNextCell(3, ids);
        ptIdx += 3;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetPolys(tris);
    tris->Delete();
    vtk_pts->Delete();

    return pd;
}

class Tetrahedron
{
  public:
    float X[4];
    float Y[4];
    float Z[4];
    float F[4];
    void PrintSelf()
    {
        for (int i = 0 ; i < 4 ; i++)
            printf("\tV%d: (%f, %f, %f) = %f\n", i, X[i], Y[i], Z[i], F[i]);
    };
};

/*
Return the case number of the cell based on the corner values and the isovalue.

cornersF = a 4 float long list that is the F values at each of the 4 corners of the cell
    Order of the values: [0] = left, [1] = right, [2] = back, [3] = top

isoval = the isovalue that's being sliced at.
*/

int IdentifyCase(float* cornersF, float isoval) {
    //categorize the cell corners as greater than or equal, or less than the isovalue.
    //more than or equal to isovalue = 1 else 0.
    //bit order: Top, Back, Right, Left
    //"bits" are converted to an int, so the bit corrisponding to "Top" by my convention would add 8 to the case num value

    int caseNum = 0;
    if (cornersF[3] >= isoval) {
        caseNum += 8;
    }

    if (cornersF[2] >= isoval) {
        caseNum += 4;
    }

    if (cornersF[1] >= isoval) {
        caseNum += 2;
    }

    if (cornersF[0] >= isoval) {
        caseNum += 1;
    }
    return caseNum;
}


/*
Interpolates along a line in 3D space to locate the point with a given isovalue
Takes the tet object which has the points as attributes

Takes the point index values based on my convention (Left = 0, Right = 1, Back = 2, Top = 3)

Takes a location to store the output point location (float[3])

Takes the isovalue
*/

void InterpEdge(Tetrahedron tet, int* points, float* pt, float isoval) {
    //t = (isovalue - F(A)) / (F(B) - F(A))
    float t;//fraction across the line
    float distX;//how much more to the right pt is from the first point
    float distY;//how much higher pt is from the first point
    float distZ;//how much further back pt is from the first point

    t = std::abs ((isoval - tet.F[points[0]]) / (tet.F[points[1]] - tet.F[points[0]]));

    //to find the total distance of line in the X/Y/Z plane, subtract tet.X/Y/Z[points[1]] from tet.X/Y/Z[points[0]]
    //then to use t to find and add the portion of that axis pt is across to tet.X/Y/Z[points[0]]. 
    //If the output of these operations is negative, then whatever. So we're moving from top to bottom instead of bottom to top.
    //My conventions pretend it's a certain case, but it often isn't.

    distX = t * (tet.X[points[1]] - tet.X[points[0]]);
    distY = t * (tet.Y[points[1]] - tet.Y[points[0]]);
    distZ = t * (tet.Z[points[1]] - tet.Z[points[0]]);

    pt[0] = tet.X[points[0]] + distX;
    pt[1] = tet.Y[points[0]] + distY;
    pt[2] = tet.Z[points[0]] + distZ;
}


void IsosurfaceTet(Tetrahedron &tet, TriangleList &tl, float isoval)
{
    /*
    Passed a tet that consists of 4 points, a triangle list to add to, and the isoval.
    
    to add a triangle to the list:   tl.addTriangle(float X1, float Y1, float Z1, 
                                                    float X2, float Y2, float Z2, 
                                                    float X3, float Y3, float Z3)

    The passed in tet consists of 4 values, X[4], Y[4], Z[4], F[4] where it contains the X/Y/Z/F values for a point when using the same index on each list.

    Convention:
    Assume the tet is placed on the table with one side directly facing you.
    The left point is 0
    The right point is 1
    The back point is 2
    The top point is 3.

    The bottom edge closeset to you is 0
    The back left edge is 1
    The back right edge is 2
    The upper left edge is 3
    the upper right edge is 4
    the upper back edge is 5

    */

    //Number of triangles in each case
    int numTris[16] = {0, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 0};

    //for every edge, what points are it connected to
    int lineEnds[6][2] = {
        {0,1},
        {0,2},
        {1,2},
        {0,3},
        {1,3},
        {2,3}
        };

    //Numbers of the edges that are touched in each case, broken into 2 packs of 3 that represent the edges touched of the triangles to render.
    //EX: If edges 1,2,3,4 have points on them, then it's lup entry would be {1,2,3,2,3,4} as this will render the full 4 point plane in the tet. 
    int lup[16][6]{
        {-1,-1,-1,-1,-1,-1},//0
        {0,1,3,-1,-1,-1},//1
        {0,2,4,-1,-1,-1},//2
        {1,2,3,2,3,4},//3
        {1,2,5,-1,-1,-1},//4
        {0,2,3,2,3,5},//5
        {0,1,4,1,4,5},//6
        {3,4,5,-1,-1,-1},//7
        //below is mirror order of above
        {3,4,5,-1,-1,-1},//8
        {0,1,4,1,4,5},//9
        {0,2,3,2,3,5},//10
        {1,2,5,-1,-1,-1},//11
        {1,2,3,2,3,4},//12
        {0,2,4,-1,-1,-1},//13
        {0,1,3,-1,-1,-1},//14
        {-1,-1,-1,-1,-1,-1} //15
        };

    //find the case of the tet: Pass in the F values list of the tet, 
    //where F[0] is the left, F[1] is the right, F[2] is the back, and F[3] is the top
    int tetCase = IdentifyCase(tet.F, isoval);

    int i;
    for (i = 0; i < numTris[tetCase]; i++) {
        //for every triangle, interpolate to the spot along the 3 edges of the triangle that the isovalues is at.
        float pt1[3];
        float pt2[3];
        float pt3[3];

        //lineEnds takes the edge number and returns the points the edge ends on.
        //lup takes the case num and then and index, and returns the edge number of that point of the triangle.

        InterpEdge(tet, lineEnds[lup[tetCase][0 + (i * 3)]], pt1, isoval);
        InterpEdge(tet, lineEnds[lup[tetCase][1 + (i * 3)]], pt2, isoval);
        InterpEdge(tet, lineEnds[lup[tetCase][2 + (i * 3)]], pt3, isoval);

        tl.AddTriangle(pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], pt3[0], pt3[1], pt3[2]);
    }
}

int main()
{
    int  i, j;

    // If you want to try on one tetrahedron.
    /*
    Tetrahedron t;
    t.X[0] = 0;//right
    t.Y[0] = 0;
    t.Z[0] = 0;
    t.F[0] = 0;

    t.X[1] = 1;//left
    t.Y[1] = 0;
    t.Z[1] = 0;
    t.F[1] = 0;

    t.X[2] = 0;//back
    t.Y[2] = 1;
    t.Z[2] = 0;
    t.F[2] = 1;

    t.X[3] = 0.5;//top
    t.Y[3] = 0.5;
    t.Z[3] = 1.0;
    t.F[3] = 0;
    TriangleList tl;
    IsosurfaceTet(t, tl, 0.5);
    */
    
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6.vtk");
    rdr->Update();
    if (rdr->GetOutput() == NULL || rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Could not find input file." << endl;
        exit(EXIT_FAILURE);
    }

    vtkUnstructuredGrid *ugrid = (vtkUnstructuredGrid *) rdr->GetOutput();
    float *pts = (float *) ugrid->GetPoints()->GetVoidPointer(0);
    float *F = (float *) ugrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    vtkCellArray *cell_array = ugrid->GetCells();
    
    TriangleList tl;
    cell_array->InitTraversal();
    int ncells = cell_array->GetNumberOfCells();
    cerr << "Number of cells to tetrahedralize is " << ncells << endl;
    
    int cnt = 0;
    
    float isoval = 12.2;
    
    
    for (int i = 0 ; i < ncells ; i++)
    {
        vtkIdType npts;
        vtkIdType *ids;
        cell_array->GetNextCell(npts, ids);
        if (npts == 4)
        {
            Tetrahedron tet;
            for (int j = 0 ; j < 4 ; j++)
            {
                // This data set is in a huge bounding box.  Normalize as we go.
                tet.X[j] = (pts[3*ids[j]]+3e+7)/6e+7;
                tet.Y[j] = (pts[3*ids[j]+1]+3e+7)/6e+7;
                tet.Z[j] = (pts[3*ids[j]+2]+3e+7)/6e+7;
                tet.F[j] = F[ids[j]];
            }
            IsosurfaceTet(tet, tl, isoval);
        }
        else
        {
            cerr << "Input was non-tetrahedron!!  Ignoring..." << endl;
            cerr << "Type is " << npts << endl;
            cerr << "Cell is " << i << endl;
            cnt++;
            continue;
        }
    }
    
    vtkPolyData *pd = tl.MakePolyData();

/*
    //This can be useful for debugging
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("proj6_out.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkCleanPolyData *cpd = vtkCleanPolyData::New();
    cpd->SetInputData(pd);
    cpd->SetAbsoluteTolerance(0);
    cpd->PointMergingOn();
    cpd->Update();
    vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
    pdn->SetInputData(cpd->GetOutput());
    //pdn->SetInputData(pd);
    pdn->Update();

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pdn->GetOutput());
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0.5, 0.5, 0.5);
    ren1->GetActiveCamera()->SetPosition(0,0,-2);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(0.01, 4);
    //ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
