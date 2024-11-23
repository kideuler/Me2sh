#include "gtest/gtest.h"
#include "geometry.hpp"

void Create_circle_coords(Me2sh_Geometry &geo, int npoints) {
    double x, y;
    for (int i = 0; i < npoints; i++) {
        x = 0.5 + 0.5*cos(2*M_PI*i/npoints);
        y = 0.5 + 0.5*sin(2*M_PI*i/npoints);
        geo.points.push_back({x, y, 0.0});
    }
}

TEST(geomTest,createCircle) {
    Me2sh_Geometry geo;
    geo.addEllipse(0.5,0.5,0.2,0.2);
    gmsh::model::mesh::generate(2);
}

TEST(geomTest,createEllipse) {
    Me2sh_Geometry geo;
    geo.addEllipse(0.5,0.5,0.3,0.2);
    gmsh::model::mesh::generate(2);
}

TEST(geomTest,createRectangle) {
    Me2sh_Geometry geo;
    geo.addRectangle(0.5,0.5,0.3,0.2);
    gmsh::model::mesh::generate(2);
}

TEST(geomTest,createSpline) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addSpline(0, np);
    gmsh::model::mesh::generate(2);
}

TEST(geomTest,createBSpline) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addBSpline(0, np);
    gmsh::model::mesh::generate(2);
}

TEST(geomTest,FuseOverlapping){
    Me2sh_Geometry geo;
    geo.addEllipse(0.5,0.5,0.1,0.1);
    geo.addEllipse(0.55,0.5,0.1,0.1);
    geo.addEllipse(0.7,0.5,0.2,0.1);
    geo.FuseOverlapping();
    gmsh::model::mesh::generate(2);
}

TEST(geomTest,MakeRectangleAndCut){
    Me2sh_Geometry geo;
    geo.addEllipse(0.5,0.5,0.1,0.1);
    geo.addEllipse(0.55,0.5,0.1,0.1);
    geo.MakeRectangleAndCut(0.5,0.5,0.5,0.5);
    gmsh::model::mesh::generate(2);
}

#ifdef USE_GEO
TEST(createBezier, geomTest) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addBezier(0, np);
    gmsh::model::mesh::generate(2);
}
#endif