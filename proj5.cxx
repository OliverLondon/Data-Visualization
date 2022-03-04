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
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
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
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>


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
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
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
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
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
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
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
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
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
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
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
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}


void BoundingBoxForCell(const float* X, const float* Y, const int* dims,
    int cellId, float* bbox)
{
    //out of range checks
    if (cellId < 0) {
        bbox[0] = -100;
        bbox[1] = 100;
        bbox[2] = -100;
        bbox[3] = 100;
        return;
    }
    if (cellId >= (dims[0] - 1) * (dims[1] - 1)) {
        bbox[0] = -100;
        bbox[1] = 100;
        bbox[2] = -100;
        bbox[3] = 100;
        return;
    }

    int column = cellId % (dims[0] - 1);
    float Float_row = cellId / (dims[0] - 1);
    int row = static_cast<int>(Float_row); //drop the decimal
    // (colulm,row)is now the logical index of the cell, assuming 0,0 is bottom left cell 

    bbox[0] = X[column];
    bbox[1] = X[column + 1];
    bbox[2] = Y[row];
    bbox[3] = Y[row + 1];
}

/*
Return the case number of the cell based on the corner values and the isovalue.

cornersF = a 4 float long list that is the F values at each of the 4 corners of the cell
    Order of the values: [0] = bottom left, [1] = bottom right, [2] = top left, [3] = top right

isoval = the isovalue that's being sliced at.
*/

int IdentifyCase(float *cornersF, float isoval) {
    //categorize the cell corners as more or less than the isovalue.
    // more than or equal to isovalue = 1 else 0.
    int cbits[4];//bit order: TR, TL, BR, BL
    
    if (cornersF[3] >= isoval) {
        cbits[0] = 1;
    }
    else { cbits[0] = 0; }

    if (cornersF[2] >= isoval) {
        cbits[1] = 1;
    }
    else { cbits[1] = 0; }

    if (cornersF[1] >= isoval) {
        cbits[2] = 1;
    }
    else { cbits[2] = 0; }

    if (cornersF[0] >= isoval) {
        cbits[3] = 1;
    }
    else { cbits[3] = 0; }

    int caseNum = 0;
    if (cbits[0]) { caseNum += 8; }
    if (cbits[1]) { caseNum += 4; }
    if (cbits[2]) { caseNum += 2; }
    if (cbits[3]) { caseNum += 1; }
    return caseNum;
}

/*
Return the X,Y value of the interpolated point along the given edge for a given isovalue.

bbox[4] = the minX, maxX, minY, maxY values of the cell the edge is a part of

edge is the edge number (bottom = 0, right = 1, top = 2, left = 3

cornerFV[4] = the F values of the corners: 
    [0] = bottom left, [1] = bottom right, [2] = top left, [3] = top right

isovalue = the given isovalue that's being interpolated to

pt[2] = the location where the X and Y coordinates are stored/returned
*/

void InterpEdgeValue(float *corners, int edge, float *cornerFV, float isovalue, float *pt) {
    //calculate T, the % of the way that the point is from left to right, or bottom to top.
    //then, once T is found, use that to interpolate the X and Y coordinate of the point
    //t = (isovalue - F(A)) / (F(B) - F(A))
    float t;//fraction across the line
    float dist;//actual distance from the left or bottom
    //size of the top/bot or left/right edges
    float topOrBot = (corners[1] - corners[0]);
    float leftOfRight = (corners[3] - corners[2]);

    if (edge == 0) {//bottom
        pt[1] = corners[2];
        t = (isovalue - cornerFV[0]) / (cornerFV[1] - cornerFV[0]);
        dist = t * topOrBot;
        pt[0] = corners[0] + dist;
    }
    if (edge == 1) {//right
        pt[0] = corners[1];
        t = (isovalue - cornerFV[1]) / (cornerFV[3] - cornerFV[1]);
        dist = t * leftOfRight;
        pt[1] = corners[2] + dist;
    }
    if (edge == 2) {//top
        pt[1] = corners[3];
        t = (isovalue - cornerFV[2]) / (cornerFV[3] - cornerFV[2]);
        dist = t * topOrBot;
        pt[0] = corners[0] + dist;
    }
    if (edge == 3) {//left
        pt[0] = corners[0];
        t = (isovalue - cornerFV[0]) / (cornerFV[2] - cornerFV[0]);
        dist = t * leftOfRight;
        pt[1] = corners[2] + dist;
    }
}

class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj5.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!
    //make and fill numSegments 
    int numSegments[16] = { 0,1,1,1,1,1,2,1,1,2,1,1,1,1,1,0 };
    /*edge/wall numbers agreed upon for the convention
    2    -2-   3

    -3-       -1-

    0    -0-   1
    
    */
    //make and fill the lookup table for edges that the lines touch for each case
    int lup[16][4] = { 
        {-1,-1,-1,-1}, {0,3,-1,-1}, {0,1,-1,-1}, {1,3,-1,-1},
        {2,3,-1,-1}, {0,2,-1,-1}, {0,3,1,2}, {1,2,-1,-1},
        {1,2,-1,-1}, {0,1,2,3}, {0,2,-1,-1}, {2,3,-1,-1},
        {1,3,-1,-1}, {0,1,-1,-1}, {0,3,-1,-1}, {-1,-1,-1,-1} };
    
    float isoval = 3.2;
    /*
    X and Y are axis lists, dims[0 or 1]  is the size of them
    F is the scalar values at each vertex,  size of F = dims[0] * dims[1] so to get the correct index:
    logicalY*dims[0] + logicalX
    */

    int nCells = GetNumberOfCells(dims);

    for (i = 0; i < nCells; i++) {
    //for (x = 0; x < 1; x++){
        //for every cell, classify it into the case x
        //then use the case to get number of line segments (nSegments) x
        //which also gives the edges that the segments end on
        //then interpolate along the edges that the segments end on.

        //get the F values at every corner
        int idx[2]; //logical X,Y coords of the cell
        GetLogicalCellIndex(idx, i, dims);
        //logical index of corners
        int BL[2] = { idx[0], idx[1] };
        int BR[2] = { idx[0]+1, idx[1] };
        int TL[2] = { idx[0], idx[1]+1 };
        int TR[2] = { idx[0]+1, idx[1]+1 };
        //F values of the corners
        float FBL = F[GetPointIndex(BL, dims)];
        float FBR = F[GetPointIndex(BR, dims)];
        float FTL = F[GetPointIndex(TL, dims)];
        float FTR = F[GetPointIndex(TR, dims)];
        float cornerFvals[4] = { FBL, FBR, FTL, FTR };

        //change the bits to an int, that is the case of the cell
        int caseNum = IdentifyCase(cornerFvals, isoval);
        int caseSegments = numSegments[caseNum];

        //similar to slide 62
        for (j = 0; j < caseSegments; j++) {
            
            int edge1 = lup[caseNum][2*j];
            int edge2 = lup[caseNum][(2*j)+1];

            //get the max and min values for X and y axis
            float bbox[4];//minX ,maxX, minY, maxY
            BoundingBoxForCell(X, Y, dims, i, bbox);

            float pt1[2];
            float pt2[2];
            //get the X,Y coords of the isovalue at the edges
            InterpEdgeValue(bbox, edge1, cornerFvals, isoval, pt1);
            InterpEdgeValue(bbox, edge2, cornerFvals, isoval, pt2);
            
            sl.AddSegment(pt1[0], pt1[1], pt2[0], pt2[1]);
        }
    }

    // -------------
    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
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

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
