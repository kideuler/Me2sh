#include "TriMesh.hpp"

int main(){
    Matrix<double> coords(1000,2);
    coords.fill_rand();
    Mesh2D msh(coords);
    msh.print();
    return 0;
}