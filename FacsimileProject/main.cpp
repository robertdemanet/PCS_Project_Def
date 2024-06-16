#include <iostream>

#include <Eigen/Geometry>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include"GeometryLibrary.hpp"

using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;


int main()
{
    string NameFile= "DFN/FR10_data.txt";
    vector<Fracture> fractures;
    fractures=readDFN(NameFile);
    vector<Trace> Traces;
    Traces=computeTraces(fractures);
    int num = fractures.size();
    string outputFileName = "Traces.txt";
    vector<vector<Support>> FractureTraces = writeResult(outputFileName,Traces,num);
    string outputFileName2 = "TracesForFracture.txt";
    bool TracesForFracture = writeTracesForFracture(outputFileName2, FractureTraces);
<<<<<<< HEAD


=======
    //PolygonalMesh mesh;
   // mesh=createMesh(fractures,mesh,FractureTraces);
>>>>>>> a9380b2c4913b4f1ab58b4bcdc6952cba8dd58ce
    return 0;
}

