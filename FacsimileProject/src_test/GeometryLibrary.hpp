#pragma once

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Eigen>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace GeometryLibrary {


// Struct per rappresentare una singola frattura
struct Fracture {
    int id;
    int numVertices;
    MatrixXd vertices;
};

// Struct per rappresentare una singola traccia
struct Trace {
    int id;
    int Fracture1ID;
    int Fracture2ID;
    Vector3d firstPoint;
    Vector3d finalPoint;
    vector<Vector3d> vertex_Inters1;
    vector<Vector3d> vertex_Inters2;
    int firstEdgeID;
    int finalEdgeID;
};

struct Support{
    int numFractures;
    int idT;
    bool Tips;
    double lenght;
    Vector3d firstPoint;
    Vector3d finalPoint;
};

// Funzione per leggere il DFN
vector<Fracture> readDFN(const string& filename);


Vector3d computeCentroid(Fracture& fracture);

bool testCircumference(Fracture& fracture1,
                       Fracture& fracture2);

vector<Trace> computeTraces(vector<Fracture>& fractures);

vector<Vector3d> TraceVertexes(Vector3d& Point1,
                               Vector3d& Point2,
                               Fracture& fracture1);

vector<double> compareAlphas(double& alpha1,double& alpha2,double& alpha3,double& alpha4);
bool compareLenght(const Support& Support1,const Support& Support2);


vector<vector<Support>> writeResult(const string& outputFilePath, // le posizioni del vector più esterno corrispondono all'id della frattura
                                    vector<Trace>& Traces,
                                    int& numFractures);              // per ogni frattura accedo al vector di struct Support e poi alla singola

bool writeTracesForFracture(const string& outputFilePath,
                            vector<vector<Support>>& FractureTraces);


// Seconda parte //////////

// CELLE 2D
struct Cell2D{
    int numIDvertices;
    int numIDedges;
    vector<unsigned int> IDs_vertices;
    vector<unsigned int> IDs_edges;
    bool status=true;
} ;

struct PolygonalMesh{
    //CELLE 0D
    vector<unsigned int> Cell0ID={};
    vector<Vector3d> Cell0DCoordinates={};


   // map<unsigned int,array< int,2>> Cell1D;//sia id del lato sia ids dei vertici
    vector<int> Cell1ID={};
    vector<Vector2i> Cell1DVertices={};

    // CELLE 2D
    vector<Cell2D> vecCell2D;

};

PolygonalMesh createMesh(vector<Fracture>& fractures,PolygonalMesh& mesh,vector<vector<Support>>& Traces);

}

















