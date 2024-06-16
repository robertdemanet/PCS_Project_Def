#ifndef __TESTUTILS_H
#define __TESTUTILS_H

#include <gtest/gtest.h>
#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include "GeometryLibrary.hpp"

using namespace std;
using namespace Eigen;

namespace GeometryLibrary{

TEST(DFNTEST,readDFNTest)
{
    // testo il numero di fratture lette
    string inputFileName = "DFN/FR3_data.txt";
    vector<Fracture> fractures = readDFN(inputFileName);
    EXPECT_EQ(fractures.size(), 3);

    // testo i dati inerenti la prima frattura:id,num di vertici e matrice di coordinate
    EXPECT_EQ(fractures[0].id, 0);
    EXPECT_EQ(fractures[0].numVertices, 4);
    MatrixXd Matrix1=MatrixXd::Zero(3, 4);
    Matrix1 << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
    EXPECT_EQ(fractures[0].vertices, Matrix1);

    // testo i dati inerenti la seconda frattura:id,num di vertici e matrice di coordinate
    EXPECT_EQ(fractures[1].id, 1);
    EXPECT_EQ(fractures[1].numVertices, 4);
    MatrixXd Matrix2=MatrixXd::Zero(3, 4);
    Matrix2 << 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01,
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
        -1.0000000000000001e-01, 2.9999999999999999e-01, 2.9999999999999999e-01, -1.0000000000000001e-01;
    EXPECT_EQ(fractures[1].vertices, Matrix2);

    // testo i dati inerenti la terza frattura:id,num di vertici e matrice di coordinate
    EXPECT_EQ(fractures[2].id, 2);
    EXPECT_EQ(fractures[2].numVertices, 4);
    MatrixXd Matrix3=MatrixXd::Zero(3, 4);
    Matrix3 << -2.3777799999999999e-01, 3.1618370000000001e-01, 3.1618370000000001e-01, -2.3777799999999999e-01,
        5.0000000000000000e-01, 5.0000000000000000e-01, 5.0000000000000000e-01, 5.0000000000000000e-01,
        -3.4444000000000002e-01, -3.4444000000000002e-01, 4.5283889999999999e-01, 4.5283889999999999e-01;
    EXPECT_EQ(fractures[2].vertices, Matrix3);
}

TEST(DFNTEST,TestComputeCentroid)
{
    // testo prendendo come frattura di esempio la prima presente nel file FR3_data
    Fracture fracture;
    fracture.numVertices = 4;
    MatrixXd Matrix=MatrixXd::Zero(3, 4);
    Matrix << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

    fracture.vertices = Matrix;
    Vector3d centroid = computeCentroid(fracture);
    Vector3d result;
    result << 0.5000000000000000e+00, 0.5000000000000000e+00,0.0000000000000000e+00;
    EXPECT_EQ(centroid, result);
}

TEST(DFNTEST,TestCircumference)
{
    // per testare prendo la rpima e la seconda frattura del file sotto
    string inputFileName = "DFN/FR3_data.txt";
    vector<Fracture> fractures = readDFN(inputFileName);
    bool var = testCircumference(fractures[0],fractures[1]);
    EXPECT_EQ(var,true);

    bool var1 = testCircumference(fractures[1],fractures[2]);
    EXPECT_EQ(var1, true);

}

TEST(DFNTEST,TestComputeTraces)
{
    Fracture fracture1;
    fracture1.numVertices = 4;
    fracture1.id = 0;
    fracture1.vertices.resize(3, 4);


    fracture1.vertices << 0, 1, 1, 0,
        0, 0, 1, 1,
        0, 0, 0, 0;


    Fracture fracture2;
    fracture2.numVertices = 4;
    fracture2.id = 1;
    fracture2.vertices.resize(3,4);
    fracture2.vertices << 0, 1, 1, 0,
        0.5, 0.5, 0.5, 0.5,
        -1, -1, 1, 1;

    Vector3d trace_vertice1;
    trace_vertice1 << 0, 0.5, 0;

    Vector3d trace_vertice2;
    trace_vertice2 << 1, 0.5, 0;


    Trace trace;

    trace.id = 0;
    trace.Fracture1ID = 0;
    trace.Fracture2ID = 1;
    trace.firstPoint = trace_vertice1;
    trace.finalPoint = trace_vertice2;
    vector<Fracture> vecFractures;
    vecFractures.push_back(fracture1);
    vecFractures.push_back(fracture2);

    vector<Trace> a;
    a=computeTraces(vecFractures);



    EXPECT_EQ(a[0].firstPoint, trace_vertice1);
    EXPECT_EQ(a[0].finalPoint, trace_vertice2);

}

TEST(DFNTEST,TestWriteResult)
{
    Vector3d trace_vertice1;
    trace_vertice1 << 0, 0.5, 0;

    Vector3d trace_vertice2;
    trace_vertice2 << 1, 0.5, 0;
    vector<Vector3d> vertex_Inters1;
    vertex_Inters1.push_back(trace_vertice2);
    vertex_Inters1.push_back(trace_vertice1);
    vector<Vector3d> vertex_Inters2;
    vertex_Inters2.push_back(trace_vertice2);
    vertex_Inters2.push_back(trace_vertice1);
    Trace trace;
    vector<Trace> vecTrace;

    trace.id = 0;
    trace.Fracture1ID = 0;
    trace.Fracture2ID = 1;
    int numFractures = 2;
    trace.firstPoint = trace_vertice1;
    trace.finalPoint = trace_vertice2;
    trace.vertex_Inters1 = vertex_Inters1;
    trace.vertex_Inters2 = vertex_Inters2;
    vecTrace.push_back(trace);
    string FileProva = "ProvaTracce.txt";
    vector<vector<Support>> r = writeResult(FileProva, vecTrace, numFractures);
    string Output = "TracesForFracture.txt";
    bool v = writeTracesForFracture(Output, r);

    EXPECT_EQ(r[0][0].idT, 0);
    EXPECT_EQ(r[0][0].Tips, false);
    EXPECT_EQ(r[0][0].lenght, 1);


    EXPECT_EQ(r[1][0].idT, 0);
    EXPECT_EQ(r[1][0].Tips, false);
    EXPECT_EQ(r[1][0].lenght, 1);

    EXPECT_EQ(v, true);
}

}

#endif
