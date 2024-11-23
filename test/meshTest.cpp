#include "gtest/gtest.h"
#include "geometry.hpp"
#include "meshing.hpp"

void Create_circle_coords(Me2sh_Geometry &geo, int npoints) {
    double x, y;
    for (int i = 0; i < npoints; i++) {
        x = 0.5 + 0.5*cos(2*M_PI*i/npoints);
        y = 0.5 + 0.5*sin(2*M_PI*i/npoints);
        geo.points.push_back({x, y, 0.0});
    }
}

TEST(meshTest,BasicTriangleMesh) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addSpline(0, np);

    Me2sh_Mesh mesh;
    mesh.generate(geo.planeTags, 6, 1, 3, 0.01);
}

TEST(meshTest,BasicQuadrilateralMesh) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addSpline(0, np);

    Me2sh_Mesh mesh;
    mesh.generate(geo.planeTags, 6, 2, 3, 0.01);
}

TEST(meshTest,DelaunayTriangleMesh) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addSpline(0, np);

    Me2sh_Mesh mesh;
    mesh.generate(geo.planeTags, 5, 1, 3, 0.01);
}

TEST(meshTest,FrontalQuadrilateralMesh) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addSpline(0, np);

    Me2sh_Mesh mesh;
    mesh.generate(geo.planeTags, 8, 2, 3, 0.01);
}

TEST(meshTest,PackParallelQuadrilateralMesh) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addSpline(0, np);

    Me2sh_Mesh mesh;
    mesh.generate(geo.planeTags, 9, 2, 3, 0.05);
}

TEST(meshTest,CrossfieldQuadrilateralMesh) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addSpline(0, np);

    Me2sh_Mesh mesh;
    mesh.generate(geo.planeTags, 11, 2, 3, 0.05);
}