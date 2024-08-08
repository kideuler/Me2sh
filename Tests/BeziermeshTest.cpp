#include "Bezier.hpp"
#include "TriMesh.hpp"

int main(){
    Matrix<double> coords(3,2);
    coords(0,0) = 0.0;
    coords(0,1) = 0.0;
    coords(1,0) = 1.0;
    coords(1,1) = 0.0;
    coords(2,0) = 0.0;
    coords(2,1) = 1.0;

    Bezier2D bezier;
    bezier.init(coords);
    
    Matrix<double> coords_msh(0,0);
    Matrix<int> segments_msh(0,0);
    bezier.create_segments(0.01, coords_msh, segments_msh);

    Mesh2D msh;
    msh.Triangulate(coords_msh,segments_msh);
    msh.Refine(0.01);
    msh.Smooth(100);

    return 0;
}