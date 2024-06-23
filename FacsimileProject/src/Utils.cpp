#include "GeometryLibrary.hpp"
#include<algorithm>
#include<cmath>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace Eigen;

namespace GeometryLibrary{

//****************************************************************************************************************

vector<Trace> computeTraces (vector<struct Fracture>& fractures)
{

    vector<Trace> vecTrace;
    Trace Trace;
    int countID=0;
    for (size_t i = 0; i < fractures.size(); ++i) {
        for (size_t j = i + 1; j < fractures.size(); ++j) {

            if (!testCircumference(fractures[i], fractures[j]))
            {
                break;
            }

            // Calcolo i vettori che generano la normale per la prima frattura:

            Vector3d u1 = fractures[i].vertices.col(2) - fractures[i].vertices.col(0);
            Vector3d v1 = fractures[i].vertices.col(1) - fractures[i].vertices.col(0);
            Vector3d norm1 = (u1.cross(v1)).normalized();

            // Calcolo i vettori che generano la normale per la seconda frattura:

            Vector3d u2 = fractures[j].vertices.col(2) - fractures[j].vertices.col(0);
            Vector3d v2 = fractures[j].vertices.col(1) - fractures[j].vertices.col(0);
            Vector3d norm2 = (u2.cross(v2)).normalized();

            // Calcolo la direzione della retta di intersezione:

            Vector3d tangent = norm1.cross(norm2);

            // Definisco il sistema lineare:

            Matrix3d A;
            A << norm1.transpose(), norm2.transpose(), tangent.transpose();

            double c=tangent.dot(tangent);

            if (A.determinant() != 0 && c!=0)
            {
                Vector3d b;
                double d1 = fractures[i].vertices.col(0).dot(norm1);
                double d2 = fractures[j].vertices.col(0).dot(norm2);
                b << d1, d2, 0;

                // Risolvo il sistema lineare:

                Vector3d intersection_point = A.fullPivLu().solve(b); //P0, utilizziamo la fattorizzazione PALU perché è la più efficiente, con un costo di (n^3)/3 operazioni
                Vector3d otherPoint = intersection_point + 0.1 * tangent;
                vector<Vector3d > vertex_Inters1 = TraceVertexes(intersection_point,otherPoint,fractures[i]); //QUI HO A,B
                vector<Vector3d > vertex_Inters2 = TraceVertexes(intersection_point,otherPoint,fractures[j]); // QUI HO C,D

                if(vertex_Inters1.size() != 2 || vertex_Inters2.size() != 2)
                {
                    continue;
                }

                double alphaA;
                double alphaB;
                double alphaC;
                double alphaD;

                if(tangent[0] != 0)
                {
                    alphaA= (vertex_Inters1[0][0]-intersection_point[0])/tangent[0];
                    alphaB= (vertex_Inters1[1][0]-intersection_point[0])/tangent[0];
                    alphaC= (vertex_Inters2[0][0]-intersection_point[0])/tangent[0];
                    alphaD= (vertex_Inters2[1][0]-intersection_point[0])/tangent[0];
                }
                else if (tangent[1] != 0)
                {
                    alphaA= (vertex_Inters1[0][1]-intersection_point[1])/tangent[1];
                    alphaB= (vertex_Inters1[1][1]-intersection_point[1])/tangent[1];
                    alphaC= (vertex_Inters2[0][1]-intersection_point[1])/tangent[1];
                    alphaD= (vertex_Inters2[1][1]-intersection_point[1])/tangent[1];
                }
                else if (tangent[2] != 0)
                {
                    alphaA= (vertex_Inters1[0][2]-intersection_point[2])/tangent[2];
                    alphaB= (vertex_Inters1[1][2]-intersection_point[2])/tangent[2];
                    alphaC= (vertex_Inters2[0][2]-intersection_point[2])/tangent[2];
                    alphaD= (vertex_Inters2[1][2]-intersection_point[2])/tangent[2];
                }

                double alpha1;
                double alpha2;
                double alpha3;
                double alpha4;

                if(alphaA>alphaB)
                {
                    alpha1=alphaB;
                    alpha2=alphaA;
                }
                else
                {
                    alpha1=alphaA;
                    alpha2=alphaB;
                }


                if(alphaC>alphaD)
                {
                    alpha3=alphaD;
                    alpha4=alphaC;
                }
                else
                {
                    alpha3=alphaC;
                    alpha4=alphaD;
                }

                vector<double> d=compareAlphas(alpha1,alpha2,alpha3,alpha4);
                if(d.size()==0)
                {
                    continue;
                }           

                Trace.vertex_Inters1=vertex_Inters1;
                Trace.vertex_Inters2=vertex_Inters2;

                Trace.id=countID;
                Trace.Fracture1ID=fractures[i].id;
                Trace.Fracture2ID=fractures[j].id;
                Trace.firstPoint=intersection_point+d[0]*tangent;
                Trace.finalPoint=intersection_point+d[1]*tangent;               

                vecTrace.push_back(Trace);
                countID=countID+1;

            }
        }
    }
    return vecTrace;
}

//****************************************************************************************************************

vector<Vector3d> TraceVertexes(Vector3d& Point1,
                               Vector3d& Point2,
                               Fracture& fracture1)
{

    const double tol= 1e-9;
    vector<Vector3d> vertex_Inters; // nelle prime due posizioni ho i punti di intersezione che trovo con la prima frattura
                                    // nelle altre due ho quelli che trovo con la frattura 2

    for(int i=0;i<fracture1.numVertices;++i)
    {
        MatrixXd A;
        A.resize(3,2);
        Vector3d b = fracture1.vertices.col(i) - Point1;

        Vector3d d = fracture1.vertices.col((i+1) % fracture1.numVertices) - fracture1.vertices.col(i); //  (i + 1) % fracture1.numVertices calcola l'indice del vertice successivo in modo ciclico.

        Vector3d t= Point2 - Point1;
        A<< t,d;

        //controllo che le rette non siano parallele
        Vector3d control = t.cross(d);
        if(control.norm() > tol)
        {
            Vector2d solution = A.householderQr().solve(b);

            Vector3d intersectionPoint = Point1 + solution[0] * (Point2 - Point1);
            // basta controllare che beta stia tra 0 e 1
            if( (-solution[1] >=0 -tol && -solution[1] <= 1+tol))
            {
                vertex_Inters.push_back(intersectionPoint);
            }
        }
    }

    return vertex_Inters;

}

//****************************************************************************************************************

vector<double> compareAlphas(double& alpha1, double& alpha2, double& alpha3, double& alpha4)
{

    if(alpha1 > alpha4 || alpha3 > alpha2)
    {
        vector<double> d = {};
        return d ;
    }

    else{
        vector<double> d = {alpha1, alpha2, alpha3, alpha4};
        sort(d.begin(),d.end());
        vector<double> d_ = {d[1], d[2]};
        return d_;
    }
}

//****************************************************************************************************************

vector<vector<Support>> writeResult(const string& outputFilePath,
                                    vector<Trace>& Traces,
                                    int& numFractures)
{
    vector<vector<Support>> Return;
    Return.resize(numFractures);
    struct Support S;
    ofstream file;
    file.open(outputFilePath);

    if (file.fail())
    {
        cerr<< "file open failed"<< endl;
    }

    const string c="; ";
    file<<"#Number of Traces: "<<endl;
    file<<Traces.size()<<endl;
    file<<"#TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<<endl;

    for(size_t i = 0; i < Traces.size(); ++i)
    {
        for(int k=0;k<3;++k)
        {
            file << Traces[i].firstPoint[k] << c;
        }
        for(int h=0;h<2;++h)
        {
            file << Traces[i].finalPoint[h] << c;
        }
        file << Traces[i].finalPoint[2] << endl;

        S.idT = Traces[i].id;
        S.lenght = pow(Traces[i].finalPoint[0]-Traces[i].firstPoint[0],2)+
                   pow(Traces[i].finalPoint[1]-Traces[i].firstPoint[1],2)+
                   pow(Traces[i].finalPoint[2]-Traces[i].firstPoint[2],2);
        //con questa condizione verifico se i vertici della traccia appartengono entrambi alla frattura1

        if(Traces[i].vertex_Inters1.size() == 2){
            if( ((Traces[i].firstPoint.isApprox(Traces[i].vertex_Inters1[0])) || (Traces[i].firstPoint.isApprox(Traces[i].vertex_Inters1[1]))) &&
                ((Traces[i].finalPoint.isApprox(Traces[i].vertex_Inters1[0])) || (Traces[i].finalPoint.isApprox(Traces[i].vertex_Inters1[1])))  )
            {
                S.Tips=false;
            }
            else
            {
                S.Tips=true;
            }
        }


        S.firstPoint=Traces[i].firstPoint;
        S.finalPoint=Traces[i].finalPoint;

        Return[Traces[i].Fracture1ID].push_back(S);

        if(Traces[i].vertex_Inters2.size()==2){
            if( ((Traces[i].firstPoint.isApprox(Traces[i].vertex_Inters2[0])) || (Traces[i].firstPoint.isApprox(Traces[i].vertex_Inters2[1]))) &&
                ((Traces[i].finalPoint.isApprox(Traces[i].vertex_Inters2[0])) || (Traces[i].finalPoint.isApprox(Traces[i].vertex_Inters2[1])))  )
            {
                S.Tips=false;
            }
            else
            {
                S.Tips=true;
            }
        }

        Return[Traces[i].Fracture2ID].push_back(S);
    }

    file.close();

    //ordino i vettori per lunghezza in maniera decrescente e separatamente per passante e non passante
    for(auto& vec:Return) // auto& deduce automaticamente il tipo degli elementi in Return e vec è una referenza a ciascun elemento.
    {
        vector<struct Support> trueSupports;
        vector<struct Support> falseSupports;
        for(auto& support:vec)
        {
            if(support.Tips)
            {
                trueSupports.push_back(support);
            }
            else
            {
                falseSupports.push_back(support);
            }

        }

        // vec è una variabile di ciclo che rappresenta un riferimento a ciascun vector<Support> contenuto in Return,
        // e viene utilizzata per riorganizzare gli elementi di questi vettori in base a specifiche condizioni e criteri di ordinamento.
        sort(trueSupports.begin(),trueSupports.end(),compareLenght);
        sort(falseSupports.begin(),falseSupports.end(),compareLenght);
        vec.clear();
        vec.insert(vec.end(),falseSupports.begin(),falseSupports.end());
        vec.insert(vec.end(),trueSupports.begin(),trueSupports.end());  // Quindi per ogni frattura ho prima tutte le tracce passanti ordinate
                                                                        // per lunghezza decrescente e poi tutte quelle non passanti ordinate
                                                                        // allo stesso modo
    }

    return Return;

}
//****************************************************************************************************************
bool compareLenght(const Support& Support1,
                   const Support& Support2)
{
    return Support1.lenght>Support2.lenght;
}
//****************************************************************************************************************
bool writeTracesForFracture(const string& outputFilePath,
                            vector<vector<Support>>& FractureTraces)
{

    ofstream file;
    file.open(outputFilePath);

    if (file.fail())
    {
        cerr<< "file open failed" << endl;
        return false;
    }

    string c="; ";
    for(size_t i = 0; i < FractureTraces.size(); ++i)
    {
        file << "# FractureId; NumTraces" << endl;
        file << i << c << FractureTraces[i].size() << endl;
        for(size_t j=0;j<FractureTraces[i].size();++j)
        {
            file << "# TraceId; Tips; Length" << endl;
            file << FractureTraces[i][j].idT << c << FractureTraces[i][j].Tips << c << FractureTraces[i][j].lenght << endl;
        }
        file << endl;
    }
    return true;
}


//****************************************************************************************************************
bool testCircumference(Fracture& fracture1,
                       Fracture& fracture2)
{

    // Definisco i centroidi per ciascuna frattura:

    Vector3d centroid1 = computeCentroid(fracture1);
    Vector3d centroid2 = computeCentroid(fracture2);

    // Inizializzo i raggi delle due circonferenze a 0:

    double radius1 = 0.00;
    double radius2 = 0.00;

    // Itero sui vertici della prima frattura, aggiornando ad ogni iterata
    // il valore massimo assunto dalla distanza del centroide da un vertice

    for (int v = 0; v < fracture1.numVertices; ++v)
    {
        Vector3d vertex1(fracture1.vertices(0,v), fracture1.vertices(1,v), fracture1.vertices(2,v));
        Vector3d diff1 = centroid1 - vertex1;
        double distance1 = diff1.norm();
        radius1 = std::max(radius1, distance1);
    }

    // Itero sui vertici della prima frattura, aggiornando ad ogni iterata
    // il valore massimo assunto dalla distanza del centroide da un vertice

    for (int v = 0; v < fracture2.numVertices; ++v)
    {
        Vector3d vertex2(fracture2.vertices(0,v), fracture2.vertices(1,v), fracture2.vertices(2,v));
        Vector3d diff2 = centroid2 - vertex2;
        double distance2 = diff2.norm();
        radius2 = std::max(radius2, distance2);
    }

    // Calcolo la distanza nella spazio 3D dei due centroidi:

    Vector3d diff_centroids = centroid1 - centroid2;
    double distance_centroids = diff_centroids.norm();

    // E verifico se essa è maggiore della somma dei due raggi:

    if(distance_centroids > (radius1 + radius2))
    {
        return false;
    }
    else
    {
        return true;
    }
}

//****************************************************************************************************************
Vector3d computeCentroid(Fracture& fracture)
{

    // Definisco il vettore (punto 3D) centroide:

    Vector3d centroid;

    // Inizializzo le somme delle coordinate x/y/z a zero:

    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;

    // Itero sul numero di vertici sommando ogni volta su ciascuna coordinata:

    int n_vert = fracture.numVertices;
    for (int j = 0; j < n_vert; ++j)
    {
        sum_x = sum_x + fracture.vertices(0,j);
        sum_y = sum_y + fracture.vertices(1,j);
        sum_z = sum_z + fracture.vertices(2,j);
    }

    // Calcolo i valori medi assunti da ciascuna coordinata e li inserisco nel vettore centroide:

    double median_x = sum_x/n_vert;
    double median_y = sum_y/n_vert;
    double median_z = sum_z/n_vert;

    centroid << median_x, median_y, median_z;
    return centroid;
}

//****************************************************************************************************************
vector<Fracture> readDFN(const string &filename) {
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error: it is impossible to open the file" << endl;
        exit(1);
    }

    string line;
    getline(file,line); // header of numFractures

    getline(file,line); // line of numFractures
    int numFractures;
    istringstream convertNF(line);
    convertNF >> numFractures;

    vector<Fracture> fractures;
    for (int i = 0; i < numFractures; ++i)
    {
        Fracture fracture;

        getline(file,line); // header of fractureID and numVertices
        getline(file,line); // line of fractureID and numVertices
        char b;
        istringstream convertIDandVERT(line);
        convertIDandVERT >> fracture.id >> b >> fracture.numVertices;
        fracture.vertices.resize(3, fracture.numVertices);

        getline(file,line); // header of vertices
        for (int j = 0; j < 3; ++j)
        {
            getline(file,line); // line of x/y/z coordinates
            char d;
            istringstream convertCOORD(line);
            for (int k = 0; k < fracture.numVertices; ++k)
            {
                if(k != fracture.numVertices - 1)
                {
                    double vertice;
                    convertCOORD >> vertice;
                    fracture.vertices(j,k) = vertice;
                    convertCOORD >> d;
                }
                else
                {
                    double vertice;
                    convertCOORD >> vertice;
                    fracture.vertices(j,k) = vertice;
                }
            }
        }
        fractures.push_back(fracture);
    }

    file.close();
    return fractures;
}
//****************************************************************************************************************

PolygonalMesh createMesh(vector<Fracture>& fractures,PolygonalMesh& mesh,vector<vector<Support>>& Traces){
    double tol = 1e-9;
    for(size_t i = 0; i < fractures.size(); ++i)
    {
        int numVertices = fractures[i].numVertices;
        for(int j = 0; j < numVertices; ++j)
        {
            mesh.Cell0ID.push_back(j);
            mesh.Cell0DCoordinates.push_back(fractures[i].vertices.col(j));

            mesh.Cell1ID.push_back(j);
            Vector2i Vertices={j, (j+1) % numVertices};
            mesh.Cell1DVertices.push_back(Vertices);
        }

        Cell2D cell2D;
        cell2D.numIDvertices= mesh.Cell1DVertices.size();
        cell2D.numIDedges=mesh.Cell1ID.size();

        for(size_t r = 0; r < mesh.Cell1DVertices.size(); ++r)
        {
            cell2D.IDs_vertices.push_back(mesh.Cell1DVertices[r][0]);

        }
        cell2D.IDs_edges.insert(cell2D.IDs_edges.end(),mesh.Cell1ID.begin(),mesh.Cell1ID.end());
        mesh.vecCell2D.push_back(cell2D);

        size_t dim=1;
        for(size_t k = 0; k < Traces[fractures[i].id].size(); ++k)
        {
            for(size_t h = 0; h < dim; ++h)
            {
                Cell2D firstCell2D;
                Cell2D secondCell2D;
                Vector3d Vertex1;
                Vector3d Vertex2;



                if(Traces[fractures[i].id][k].finalPoint[0] == Traces[fractures[i].id][k].firstPoint[0] &&
                   Traces[fractures[i].id][k].finalPoint[1] < Traces[fractures[i].id][k].firstPoint[1] )  // il punto con la coordinata y più bassa (quando le coordinate x sono uguali) viene assegnato a Vertex1
                {
                    Vertex1=Traces[fractures[i].id][k].finalPoint;
                    Vertex2=Traces[fractures[i].id][k].firstPoint;
                }
                else
                {
                    Vertex1=Traces[fractures[i].id][k].firstPoint;
                    Vertex2=Traces[fractures[i].id][k].finalPoint;
                }

                if(mesh.vecCell2D[h].status==true)
                {
                    vector<unsigned int> edges = mesh.vecCell2D[h].IDs_edges;

                    unsigned int id_Vertex1 = numVertices;
                    unsigned int id_Vertex2 = numVertices+1;



                    if(Traces[fractures[i].id][k].Tips==false) //CONTROLLARE SE PRENDE FALSE O 1
                    {

                        int index1;
                        int index2;
                        for(int z=0;z<mesh.vecCell2D[h].numIDedges;++z)
                        {
                            Vector3d Coordinate1 = mesh.Cell0DCoordinates[mesh.vecCell2D[h].IDs_vertices[z]];
                            Vector3d Coordinate2 = mesh.Cell0DCoordinates[mesh.vecCell2D[h].IDs_vertices[(z+1) % mesh.vecCell2D[h].IDs_vertices.size()]];

                            Vector3d AB = Coordinate2-Coordinate1;
                            Vector3d AP = Vertex1-Coordinate1;
                            Vector3d AP2 = Vertex2-Coordinate1;
                            Vector3d t = AB.cross(AP);
                            Vector3d t2 = AB.cross(AP2);
                            // Il prodotto vettoriale di due vettori sarà nullo (un vettore di tutti zeri)
                            // se e solo se i due vettori sono paralleli, cioè se sono collineari.
                            if(abs(t[0]) < tol && abs(t[1]) < tol && abs(t[2]) < tol)
                            {
                                index1 = z;
                            }


                            if(abs(t2[0]) < tol && abs(t2[1]) < tol && abs(t2[2]) < tol)
                            {
                                index2 = z;
                            }
                        }




                        Vector2i vec={id_Vertex2,id_Vertex1};
                        Vector2i val={id_Vertex1,mesh.Cell1DVertices[edges[index1]][1]};
                        Vector2i val2={id_Vertex2,mesh.Cell1DVertices[edges[index2]][1]};

                        // Aggiorno le Celle0D
                        mesh.Cell0ID.push_back(id_Vertex1);
                        mesh.Cell0ID.push_back(id_Vertex2);
                        mesh.Cell0DCoordinates.push_back(Vertex1);
                        mesh.Cell0DCoordinates.push_back(Vertex2);

                        // INDEX 1
                        //AGGIORNO LE CELLE 1D
                        mesh.Cell1DVertices[edges[index1]][1]=id_Vertex1;
                        mesh.Cell1DVertices.insert(mesh.Cell1DVertices.begin()+ edges[index1]+1,val);
                        mesh.Cell1ID.push_back(numVertices);

                        //AGGIORNO LE CELLE 2D
                        mesh.vecCell2D[h].IDs_vertices.insert(mesh.vecCell2D[h].IDs_vertices.begin()+index1+1,id_Vertex1);
                        mesh.vecCell2D[h].numIDvertices=mesh.vecCell2D[h].numIDvertices+1;
                        mesh.vecCell2D[h].IDs_edges.push_back(numVertices);
                        mesh.vecCell2D[h].numIDedges=mesh.vecCell2D[h].numIDedges+1;


                        //INDEX 2
                        //AGGIORNO LE CELLE 1D
                        mesh.Cell1DVertices[edges[index2]+1][1]=id_Vertex2;
                        mesh.Cell1DVertices.insert(mesh.Cell1DVertices.begin()+ edges[index2]+2,val2);
                        mesh.Cell1DVertices.push_back(vec);
                        mesh.Cell1ID.push_back(numVertices+1);

                        //AGGIORNO LE CELLE 2D
                        mesh.vecCell2D[h].IDs_vertices.insert(mesh.vecCell2D[h].IDs_vertices.begin()+index2+2,id_Vertex2);
                        mesh.vecCell2D[h].numIDvertices=mesh.vecCell2D[h].numIDvertices+1;
                        mesh.vecCell2D[h].IDs_edges.push_back(numVertices+1);
                        mesh.vecCell2D[h].numIDedges=mesh.vecCell2D[h].numIDedges+1;
                    }

                    //tracce non passanti
                    else
                    {
                        int index1;
                        int index2;
                        vector<Vector3d> vertex_Inters;
                        vector<int> indexes;
                        Vector3d c = Vertex2 - Vertex1;

                        for(int z = 0; z < mesh.vecCell2D[h].numIDedges; ++z)
                        {
                            Vector3d Coordinate1 = mesh.Cell0DCoordinates[mesh.vecCell2D[h].IDs_vertices[z]];
                            Vector3d Coordinate2 = mesh.Cell0DCoordinates[mesh.vecCell2D[h].IDs_vertices[(z+1) % mesh.vecCell2D[h].IDs_vertices.size()]];
                            MatrixXd A(3,2);
                            Vector3d b = Coordinate1 - Vertex1;
                            Vector3d d = Coordinate2 - Coordinate1;

                            A<< c,d;

                            //controllo che le rette non siano parallele
                            Vector3d control = c.cross(d);

                            if(control.norm() > tol)
                            {
                                Vector2d solution = A.householderQr().solve(b);
                                Vector3d intersectionPoint = Vertex1 + solution[0] * c;

                                if( (-solution[1] >= 0 - tol && -solution[1] <= 1 + tol))
                                {
                                    vertex_Inters.push_back(intersectionPoint);
                                    indexes.push_back(z);
                                }
                            }
                        }


                        Vector3d vertex1 = vertex_Inters[0];
                        Vector3d vertex2 = vertex_Inters[1];

                        index1 = indexes[0];
                        index2 = indexes[1];
                        Vector2i vec = {id_Vertex2,id_Vertex1};
                        Vector2i val = {id_Vertex1,mesh.Cell1DVertices[edges[index1]][1]};
                        Vector2i val2 = {id_Vertex2,mesh.Cell1DVertices[edges[index2]][1]};
                        mesh.Cell0ID.push_back(id_Vertex1);
                        mesh.Cell0ID.push_back(id_Vertex2);
                        mesh.Cell0DCoordinates.push_back(vertex1);
                        mesh.Cell0DCoordinates.push_back(vertex2);

                        // INDEX 1
                        //AGGIORNO LE CELLE 1D
                        mesh.Cell1DVertices[edges[index1]][1] = id_Vertex1;
                        mesh.Cell1DVertices.insert(mesh.Cell1DVertices.begin() + edges[index1] + 1, val);
                        mesh.Cell1ID.push_back(numVertices);

                        //AGGIORNO LE CELLE 2D
                        mesh.vecCell2D[h].IDs_vertices.insert(mesh.vecCell2D[h].IDs_vertices.begin() + index1 + 1, id_Vertex1);
                        mesh.vecCell2D[h].numIDvertices = mesh.vecCell2D[h].numIDvertices + 1;
                        mesh.vecCell2D[h].IDs_edges.push_back(numVertices);
                        mesh.vecCell2D[h].numIDedges = mesh.vecCell2D[h].numIDedges + 1;

                        //INDEX 2
                        //AGGIORNO LE CELLE 1D
                        mesh.Cell1DVertices[edges[index2]+1][1] = id_Vertex2;
                        mesh.Cell1DVertices.insert(mesh.Cell1DVertices.begin() + edges[index2] + 2, val2);
                        mesh.Cell1DVertices.push_back(vec);
                        mesh.Cell1ID.push_back(numVertices+1);

                        //AGGIORNO LE CELLE 2D
                        mesh.vecCell2D[h].IDs_vertices.insert(mesh.vecCell2D[h].IDs_vertices.begin() + index2 + 2, id_Vertex2);
                        mesh.vecCell2D[h].numIDvertices=mesh.vecCell2D[h].numIDvertices+1;
                        mesh.vecCell2D[h].IDs_edges.push_back(numVertices+1);
                        mesh.vecCell2D[h].numIDedges=mesh.vecCell2D[h].numIDedges+1;
                    }

                    auto initV1 =find(mesh.vecCell2D[h].IDs_vertices.begin(),mesh.vecCell2D[h].IDs_vertices.end(),id_Vertex1);
                    // find è utilizzata per cercare il primo elemento in un intervallo che sia uguale a un valore specificato
                    // auto permette al compilatore di dedurre automaticamente il tipo della variabile initV1 in base al tipo di ritorno di find.
                    size_t index=distance(mesh.vecCell2D[h].IDs_vertices.begin(),initV1);

                    // itera attraverso i vertici di mesh.vecCell2D[h].IDs_vertices finché non trova id_Vertex2.
                    while(mesh.vecCell2D[h].IDs_vertices[index % mesh.vecCell2D[h].IDs_vertices.size()]!=id_Vertex2)
                    {
                        firstCell2D.IDs_vertices.push_back(mesh.vecCell2D[h].IDs_vertices[index % mesh.vecCell2D[h].IDs_vertices.size()]);
                        firstCell2D.IDs_edges.push_back(index % mesh.vecCell2D[h].IDs_vertices.size());
                        index = index + 1;
                    }
                    // Dopo che il primo ciclo si è fermato su id_Vertex2, il secondo ciclo continua dal punto in cui
                    // si è fermato il primo ciclo e itera fino a trovare id_Vertex1
                    while(mesh.vecCell2D[h].IDs_vertices[index % mesh.vecCell2D[h].IDs_vertices.size()]!=id_Vertex1)
                    {
                        secondCell2D.IDs_vertices.push_back(mesh.vecCell2D[h].IDs_vertices[index % mesh.vecCell2D[h].IDs_vertices.size()]);
                        secondCell2D.IDs_edges.push_back(index % mesh.vecCell2D[h].IDs_vertices.size());
                        index = index + 1;
                    }

                    firstCell2D.numIDvertices=firstCell2D.IDs_vertices.size()+1;
                    firstCell2D.numIDedges=firstCell2D.IDs_edges.size()+1;
                    firstCell2D.IDs_edges.push_back(mesh.Cell1ID.size());
                    firstCell2D.IDs_vertices.push_back(id_Vertex2);


                    secondCell2D.numIDvertices=secondCell2D.IDs_vertices.size()+1;
                    secondCell2D.numIDedges=secondCell2D.IDs_edges.size()+1;
                    secondCell2D.IDs_edges.push_back(mesh.Cell1ID.size());
                    secondCell2D.IDs_vertices.push_back(id_Vertex1);

                    mesh.vecCell2D[h].status=false;
                    mesh.vecCell2D.push_back(firstCell2D);
                    mesh.vecCell2D.push_back(secondCell2D);
                    numVertices=numVertices+2;
                }
            }
            dim=mesh.vecCell2D.size();
        }
    }
    return mesh;
}

}
