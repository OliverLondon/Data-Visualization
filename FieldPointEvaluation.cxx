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
#include <vtkFloatArray.h>


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
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)*idx[0];
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

// ****************************************************************************
//  Function: BoundingBoxForCell
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     cellId: a cellIndex (I.e., between 0 and GetNumberOfCells(dims))
//     bbox (output): the bounding box of cellId.  Format should be
//                     bbox[0]: the minimum X value in cellId.
//                     bbox[1]: the maximum X value in cellId.
//                     bbox[2]: the minimum Y value in cellId.
//                     bbox[3]: the maximum Y value in cellId.
//
//  Returns:  None (argument bbox is output)
//
// ****************************************************************************

void
BoundingBoxForCell(const float* X, const float* Y, const int* dims,
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
    if (cellId >= (dims[0] - 1)* (dims[1] - 1)) {
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
    // IMPLEMENT ME!
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float
EvaluateFieldAtLocation(const float* pt, const int* dims,
    const float* X, const float* Y, const float* F)
{
    if (pt[0] < X[0]) {//X too low
        return 0;
    }
    if (pt[0] > X[dims[0] - 1]) {//X too high
        return 0;
    }
    if (pt[1] < Y[0]) {//Y too low
        return 0;
    }
    if (pt[1] > Y[dims[1] - 1]) {//Y too high
        return 0;
    }

    int logicalX, logicalY, i;

    //find the logical cell index the point is in
    for (i = 1; i < dims[0]; i++) {//check X
        if (pt[0] < X[i]) {
            logicalX = i - 1;
            break;
        }
    }
    for (i = 1; i < dims[1]; i++) {//check Y
        if (pt[1] < Y[i]) {
            logicalY = i - 1;
            break;
        }
    }

    int cellId = (logicalY * (dims[0] - 1)) + logicalX; //get cellId for bounding call
    float bbox[4];
    BoundingBoxForCell(X, Y, dims, cellId, bbox);

    //fraction of distance to maxX from minX that pt[0] is at
    float distToBFromA_X = bbox[1] - bbox[0]; //right - left
    float distToPTFromA_X = pt[0] - bbox[0]; //pointX - left
    float fractionAcross_X = distToPTFromA_X / distToBFromA_X; //fraction from A to B = t
    //fraction of distance to maxY from minY that pt[1] is at
    float distToBFromA_Y = bbox[3] - bbox[2]; //top - bottom
    float distToPTFromA_Y = pt[1] - bbox[2]; //pointY - bottom
    float fractionAcross_Y = distToPTFromA_Y / distToBFromA_Y; //fraction from A to B = t

    //formula: F(A) + (t * (F(B) - F(A))) where A is bottom or left, and B is top or right depending if its X or Y axis in question.
	//									F(A)									t							F(B)									F(A)
    float interpBottomXval = F[logicalX + (logicalY * dims[0])]       + (fractionAcross_X * (F[logicalX + (logicalY * dims[0]) + 1]       - F[logicalX + (logicalY * dims[0])]));
    float interpTopXval    = F[logicalX + ((logicalY + 1) * dims[0])] + (fractionAcross_X * (F[logicalX + ((logicalY + 1) * dims[0]) + 1] - F[logicalX + ((logicalY + 1) * dims[0])]));

    //get the value of it in the y axis, which will be the final value
    float interpYval = interpBottomXval + (fractionAcross_Y * (interpTopXval - interpBottomXval));

    return interpYval; // IMPLEMENT ME!!
}

// ****************************************************************************
//  Function: CountNumberOfStraddingCells
//
//  Arguments:
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//  Returns:  the number of cells that straddle 0, i.e., the number of cells
//            that contains points who have F>0 and also have points with F<0.
//
// ****************************************************************************

int
CountNumberOfStraddlingCells(const float* X, const float* Y, const int* dims,
    const float* F)
{
    int straddleCount = 0;
    int cells = GetNumberOfCells(dims);
    int cellpos_logical[2];
    int cbox[4]; //point indexes for the corners
    float vbox[4]; //F values in the corners
    int logical_x, logical_y, i, j, neg_flag, pos_flag;
    //logical indexes for corners of cell
    int Corner_BL[2];
    int Corner_BR[2];
    int Corner_TL[2];
    int Corner_TR[2];

    //check every cell to see if the scalar values at the corners contain both a positive and negative number.
    for (i = 0; i < cells; i++) {
        //reset/clear values
        cellpos_logical[0] = 0;
        cellpos_logical[1] = 0;
        logical_x = 0;
        logical_y = 0;

        GetLogicalCellIndex(cellpos_logical, i, dims);//get logical id of the current cell
        logical_x = cellpos_logical[0];
        logical_y = cellpos_logical[1];

        //get corner logical indexes
        Corner_BL[0] = logical_x;
        Corner_BL[1] = logical_y;
        Corner_BR[0] = logical_x + 1;
        Corner_BR[1] = logical_y;
        Corner_TL[0] = logical_x;
        Corner_TL[1] = logical_y + 1;
        Corner_TR[0] = logical_x + 1;
        Corner_TR[1] = logical_y + 1;

        //get non-logical indexes for corners
        cbox[0] = GetPointIndex(Corner_BL, dims);
        cbox[1] = GetPointIndex(Corner_BR, dims);
        cbox[2] = GetPointIndex(Corner_TL, dims);
        cbox[3] = GetPointIndex(Corner_TR, dims);

        //use non-logical indexes to get F values of the corners
        vbox[0] = F[cbox[0]]; //bottom_left
        vbox[1] = F[cbox[1]]; //bottom_right
        vbox[2] = F[cbox[2]]; //top_left
        vbox[3] = F[cbox[3]]; //top_right

        neg_flag = 0;
        pos_flag = 0;
        for (j = 0; j < 4; j++) {
            if (vbox[j] < 0) {//if there's a negative
                neg_flag = 1;
                continue;
            }
            if (vbox[j] > 0) {//if there's a positive
                pos_flag = 1;
                continue;
            }
        }
        if (pos_flag == 1) {//if there's both
            if (neg_flag == 1) {
                straddleCount += 1;//increment counter
            }
        }
    }

    return straddleCount;
    // IMPLEMENT ME!
}

int main()
{
    int  i;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj2_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int numCells = CountNumberOfStraddlingCells(X, Y, dims, F);
    cerr << "The number of cells straddling zero is " << numCells << endl;

    float bbox[4];
    const int ncells = 5;
    int cellIds[ncells] = { 0, 50, 678, 1000, 1200 };
    for (i = 0 ; i < ncells ; i++)
    {
        BoundingBoxForCell(X, Y, dims, cellIds[i], bbox);
        cerr << "The bounding box for cell " << cellIds[i] << " is " 
             << bbox[0] << "->" << bbox[1] << ", " << bbox[2] << "->" << bbox[3]
             << endl;
    }

    const int npts = 10;
    float pt[npts][3] = 
         {
            {1.01119, 0.122062, 0},
            {0.862376, 1.33839, 0},
            {0.155026, 0.126123, 0},
            {0.69736, 0.0653565, 0},
            {0.2, 0.274117, 0},
            {0.893699, 1.04111, 0},
            {0.608791, -0.0533753, 0},
            {1.00543, 0.138024, 0},
            {0.384128, -0.0768977, 0},
            {0.666757, 0.60259, 0},
         };

    

    for (i = 0 ; i < npts ; i++)
    {
        float f = EvaluateFieldAtLocation(pt[i], dims, X, Y, F);
        cerr << "Evaluated field at (" << pt[i][0] <<"," << pt[i][1] << ") as "
             << f << endl;
    }
    
   
    cerr << "Infinite loop here, else Windows people may have the terminal "
         << "disappear before they see the output."
         << " Remove these lines if they annoy you." << endl;
    cerr << "(press Ctrl-C to exit program)" << endl;
    while (1) ; 
}




