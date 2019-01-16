/*=========================================================================
Program:   Visualization Toolkit
Module:    GPURenderDemo.cxx
Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
// VTK includes
#include "vtkBoxWidget.h"
#include "vtkCamera.h"
#include "vtkCommand.h"
#include "vtkColorTransferFunction.h"
#include "vtkDICOMImageReader.h"
#include "vtkImageData.h"
#include "vtkImageResample.h"
#include "vtkMetaImageReader.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPlanes.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkXMLImageDataReader.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkInteractorStyle.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkSmartPointer.h"
#include "vtkIntArray.h"
#include "vtkPoints.h"
#include "vtkMetaImageWriter.h"
#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkSphereSource.h"
#include "vtkActor.h"
#include "vtkLineSource.h"
// stdlib include
#include <stdlib.h>
#include <vector>
#include <list>
#include <array>

double* CalculateCylinderEdge(double* placedPoint, vtkImageData* input);
double CalculateDistance(double*, double*);
int CalculateMinDistancePoint(vtkImageData* input, std::vector<int> points, double* point, int selectedPoint, double *candidatePoint);

int main(int argc, char *argv[])
{
    // Read the data
    vtkImageData *input = nullptr;
    vtkPointData *pDat = nullptr;
    vtkDataArray *dArr = nullptr;
    double pointLabelValue = 0.0;
    vtkMetaImageReader *metaReader = vtkMetaImageReader::New();
    metaReader->SetFileName("C:\\users\\danie\\Documents\\CTLabelMaps\\Tube1.mha");
    metaReader->Update();
    input = metaReader->GetOutput();

    pDat = input->GetPointData();
    dArr = pDat->GetScalars();
    double *a = new double[3];
    double count = 0;
    double *placedPoint1 = new double[3]{ 33.405, -17.629, 5.402 };

    double *tube1Origin1 = CalculateCylinderEdge(placedPoint1, input);


    double *placedPoint2 = new double[3]{ -32.541, 13.912, 3.072 };
    double *candidatePoint = new double[3];
    int selectedPoint = 0;
    input = metaReader->GetOutput();
    metaReader->Update();



}

// 
double* CalculateCylinderEdge(double* placedPoint, vtkImageData* input)
{

    double pointLabelValue = 0.0;
    vtkMetaImageReader *metaReader = vtkMetaImageReader::New();
    vtkPointData *pDat = nullptr;
    vtkDataArray *dArr = nullptr;
    int numPoints = input->GetNumberOfPoints();

    dArr = input->GetPointData()->GetScalars();
    //dArr = pDat->GetScalars();
    double *a = new double[3];
    double count = 0;
    double *candidatePoint = new double[3];
    int selectedPoint = 0;

    std::vector<int> points;

    // Calculate value for initial minimum distance
    double minDistance = CalculateDistance(input->GetPoint(0), placedPoint);
    double dist = 0;

    // Find the closest tube point from the inside the tube
    for (int i = 1; i < input->GetNumberOfPoints(); i++)
    {
        pointLabelValue = dArr->GetVariantValue(i).ToDouble();
        //dArr->GetVariantValue(i);

        if (pointLabelValue == 1.0)
        {
            points.push_back(i);
            candidatePoint = input->GetPoint(i);
            dist = CalculateDistance(candidatePoint, placedPoint);
            if (dist < minDistance)
            {
                minDistance = dist;
                selectedPoint = i;
            }

        }

    }

    printf("%f %f %f\n", input->GetPoint(selectedPoint)[0], input->GetPoint(selectedPoint)[1], input->GetPoint(selectedPoint)[2]);
    // Find Normalized vector

    double* norm = new double[3];
    double normal = 0;
    norm[0] = placedPoint[0] - input->GetPoint(selectedPoint)[0];
    norm[1] = placedPoint[1] - input->GetPoint(selectedPoint)[1];
    norm[2] = placedPoint[2] - input->GetPoint(selectedPoint)[2];

    normal = sqrt(placedPoint[0] * placedPoint[0] + placedPoint[1] * placedPoint[1] + placedPoint[2] * placedPoint[2]);
    norm[0] = norm[0] / normal;
    norm[1] = norm[1] / normal;
    norm[2] = norm[2] / normal;
    double result = 0.0;
    double min2 = CalculateDistance(input->GetPoint(0), input->GetPoint(selectedPoint));

    bool intersectionFound = false;
    int N = 0;
    int selectedPoint2 = 0;
    double *vectorPoint = new double[3];
    double *asd = new double[3];

    ofstream myfile;
    myfile.open("C:\\users\\danie\\documents\\test1.fcsv");

    for (int N = 0; N < 10; N++)
    {
        vectorPoint[0] = input->GetPoint(selectedPoint)[0] + N*(norm[0]);
        vectorPoint[1] = input->GetPoint(selectedPoint)[1] + N*(norm[1]);
        vectorPoint[2] = input->GetPoint(selectedPoint)[2] + N*(norm[2]);
        myfile << -vectorPoint[0] << "," << -vectorPoint[1] << "," << vectorPoint[2] <<"\n";

    }
    myfile.close();

    while (!intersectionFound)
    {
        // Get next point along vector
        vectorPoint[0] = input->GetPoint(selectedPoint)[0] + N*(norm[0]);
        vectorPoint[1] = input->GetPoint(selectedPoint)[1] + N*(norm[1]);
        vectorPoint[2] = input->GetPoint(selectedPoint)[2] + N*(norm[2]);

        selectedPoint2 = CalculateMinDistancePoint(input, points, vectorPoint, selectedPoint, asd);

        if (selectedPoint2 != selectedPoint)
        {
            break;
        }
        N++;
    }
  

    printf("i: %d %f %f %f\n", selectedPoint2, input->GetPoint(selectedPoint2)[0], input->GetPoint(selectedPoint2)[1], input->GetPoint(selectedPoint2)[2]);

    return input->GetPoint(selectedPoint2);
}

double CalculateDistance(double* candidatePoint, double* placedPoint)
{
    double *vector = new double[3];
    vector[0] = candidatePoint[0] - placedPoint[0];
    vector[1] = candidatePoint[1] - placedPoint[1];
    vector[2] = candidatePoint[2] - placedPoint[2];

    double dist = abs(sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]));
    delete vector;
    return dist;
}

// Takes in a list of Id's corresponding to Tube points inside an imagedata object
// and returns the id of the point which has the closest distance to the candidate point
int CalculateMinDistancePoint(vtkImageData* input, std::vector<int> points, double* point, int selectedPoint, double* candidatePoint)
{

    double minDistance = 1000;
    double dist = 0;
    int selectedPoint2 = 0;

    // Loop through points belonging to tube label
    for (int i = 0; i < points.size(); i++)
    {
        candidatePoint = input->GetPoint(points.at(i));
        dist = CalculateDistance(candidatePoint, point);
        if (dist < minDistance)
        {
            minDistance = dist;
            selectedPoint2 = points.at(i);
        }
    }
    // index belonging to vtkImageData
    return selectedPoint;
}
