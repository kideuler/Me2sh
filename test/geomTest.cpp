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

TEST(createCircle, geomTest) {
    Me2sh_Geometry geo;
    geo.addEllipse(0.5,0.5,0.2,0.2);
    gmsh::model::mesh::generate(2);
}

TEST(createEllipse, geomTest) {
    Me2sh_Geometry geo;
    geo.addEllipse(0.5,0.5,0.3,0.2);
    gmsh::model::mesh::generate(2);
}

TEST(createSpline, geomTest) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addSpline(0, np);
    gmsh::model::mesh::generate(2);
}

TEST(createBSpline, geomTest) {
    int np = 100;
    Me2sh_Geometry geo;
    Create_circle_coords(geo, np);
    geo.addBSpline(0, np);
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