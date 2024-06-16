#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolygonalLibrary {

struct PolygonalMesh
{
    unsigned int NumberCell0D = 0; // number of Cell0D
    std::vector<unsigned int> Cell0DId = {}; // vector of size 1xNumberCell0D, containing Cell0Ds' id
    std::vector<Vector2d> Cell0DCoordinates = {}; // vector of size 2xNumberCell0D, containing Cell0Ds' coordinates (x,y)
    std::map<unsigned int, list<unsigned int>> Cell0DMarkers = {}; // vector of size 1xNumberCell0D, containing Cell0D's  markers

    unsigned int NumberCell1D = 0; // number of Cell1D
    std::vector<unsigned int> Cell1DId = {}; // vector of size 1xNumberCell1D, containing Cell1Ds' id
    std::vector<Vector2i> Cell1DVertices = {}; // vector of integers of size 2xNumberCell1D, containing Cell1Ds' vertices (fromId, toId)
    std::map<unsigned int, list<unsigned int>> Cell1DMarkers = {}; // vector of size 1xNumberCell1D, containing Cell1D's properties (marker)

    unsigned int NumberCell2D = 0; // number of Cell2D
    std::vector<unsigned int> Cell2DId = {}; // vector of size 1xNumberCell2D, containing Cell2Ds' id
    std::vector<vector<unsigned int>> Cell2DVertices = {}; // vector of size 1xNumVerticesCell2D, containing Cell2Ds' vertices
    std::vector<vector<unsigned int>> Cell2DEdges = {}; // vector of size 1xNumEdgesCell2D, containing Cell2Ds' edges
    std::map<unsigned int, list<unsigned int>> Cell2DMarkers = {}; // vector of size 1xNumberCell2D, containing Cell2D's markers
};

}
