#include "TriMesh.hpp"

int TEST1(){
    // random coordinates
    Matrix<double> coords(1000,2);
    coords.fill_rand();
    Mesh2D msh;

    // triangulating cooreinates
    msh.Triangulate(coords);
    msh.Refine(0.01);
    msh.Smooth(100);
    return 0;
}

int TEST2(){
    IntMatrix segments;
    DoubleMatrix coords;
    double h = Flower(&segments, &coords, 150, false);

    Matrix<double> coords_msh(150,2);
    Matrix<int> segments_msh(150,2);

    for (int i = 0; i<150; i++){
        coords_msh(i,0) = coords.data[2*i];
        coords_msh(i,1) = coords.data[2*i+1];
        segments_msh(i,0) = segments.data[2*i];
        segments_msh(i,1) = segments.data[2*i+1];
    }

    Mesh2D msh;

    // triangulating cooreinates
    msh.Triangulate(coords_msh,segments_msh);
    msh.Refine(0.01);
    msh.Smooth(100);
    return 0;
}

int main(){
    TEST1();
    TEST2();

    return 0;
}