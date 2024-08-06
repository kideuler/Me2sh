#include "TriMesh.hpp"

Mesh2D::Mesh2D(Matrix<double> &coords){

    int rows = coords.nrows();
    int cols = coords.ncols();
    DoubleMatrix xs = DoubleMatrix_create(rows, cols);

    for (int i = 0; i<rows; ++i){
        for (int j = 0; j<cols; ++j){
            xs.data[i*cols + j] = coords(i,j);
        }
    }

    Mesh msh = GeoMesh_Delaunay(&xs,2);
    Mesh2vtk(&msh);
}