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


int main(){
    TEST1();

    return 0;
}