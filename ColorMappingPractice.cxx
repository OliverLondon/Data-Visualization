using namespace std;
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <algorithm>


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

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F)
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
    float interpBottomXval = F[logicalX + (logicalY * dims[0])] + (fractionAcross_X * (F[logicalX + (logicalY * dims[0]) + 1] - F[logicalX + (logicalY * dims[0])]));
    float interpTopXval = F[logicalX + ((logicalY + 1) * dims[0])] + (fractionAcross_X * (F[logicalX + ((logicalY + 1) * dims[0]) + 1] - F[logicalX + ((logicalY + 1) * dims[0])]));

    //get the value of it in the y axis, which will be the final value
    float interpYval = interpBottomXval + (fractionAcross_Y * (interpTopXval - interpBottomXval));

    return interpYval;
}


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128) 
//        F=1: (255,255,255) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************


void
ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    RGB[0] = 255 * F;
    RGB[1] = 255 * F;
    RGB[2] = 128 + (127 * F);
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F=0: (0,0,128) 
//        F=0.5: (255,255,255) 
//        F=1: (128, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    float reNormF;

    if (F >= 0.5) {
        reNormF = (F - 0.5) / (1.0 - 0.5);//fraction across 0.5 -> 1
        RGB[0] = 255 - (reNormF * 127);
        RGB[1] = 255 - (reNormF * 255);
        RGB[2] = 255 - (reNormF * 255);
    }
    else {
        //(F - 0) / (0.5 - 0)
        reNormF = F * 2; //fraction across 0 -> 0.5
        RGB[0] = reNormF * 255;
        RGB[1] = reNormF * 255;
        RGB[2] = 128 + (reNormF * 127);
    }
}

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyHSVColorMap(float F, unsigned char* RGB)
{
    //from color lecture
    //hsvToRGB(float hue, float saturation, float value)
    float v = 255;
    float hue = 360.0 * F;
    hue /= 60.f;
    // sector 0 to 5
    int i = floor(hue);
    float f = hue - i;
    // factorial part of h
    float p = 0;
    float q = (1 - (1 * f)) * 255;
    float t = (1 - (1 * (1 - f))) * 255;
    switch (i) {
        case 0:
            RGB[0] = v;
            RGB[1] = t;
            RGB[2] = p;
            break;
        case 1:
            RGB[0] = q;
            RGB[1] = v;
            RGB[2] = p;
            break;
        case 2:
            RGB[0] = p;
            RGB[1] = v;
            RGB[2] = t;
            break;
        case 3:
            RGB[0] = p;
            RGB[1] = q;
            RGB[2] = v;
            break;
        case 4:
            RGB[0] = t;
            RGB[1] = p;
            RGB[2] = v;
            break;
        case 5:
            RGB[0] = v;
            RGB[1] = p;
            RGB[2] = q;
            break;
    }
}

int main()
{
    int  i, j;

    vtkDataSetReader* rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid* rgrid = (vtkRectilinearGrid*)rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float* X = (float*)rgrid->GetXCoordinates()->GetVoidPointer(0);
    float* Y = (float*)rgrid->GetYCoordinates()->GetVoidPointer(0);
    float* F = (float*)rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    int nx = 500;
    int ny = 500;

    vtkImageData* images[3];
    unsigned char* buffer[3];
    for (i = 0; i < 3; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char*)images[i]->GetScalarPointer(0, 0, 0);
    }

    for (i = 0; i < 3 * nx * ny; i++)
        for (j = 0; j < 3; j++)
            buffer[j][i] = 0;

    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            float stepval = 18.0 / 499.0;
            pt[0] = (stepval * i) - 9.0;
            pt[1] = (stepval * j) - 9.0;
            
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);

            float normalizedF = (f - 1.2) / (5.02 - 1.2);

            // I TAKE OVER HERE
            int offset = 3 * (j * nx + i);
            ApplyBlueHotColorMap(normalizedF, buffer[0] + offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1] + offset);
            ApplyHSVColorMap(normalizedF, buffer[2] + offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
