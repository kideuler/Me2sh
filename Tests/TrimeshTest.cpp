#include "TriMesh.hpp"

int main(){
    Matrix<double> coords(100,2);
    coords.fill_rand();
    Mesh2D msh(coords);
    return 0;
}