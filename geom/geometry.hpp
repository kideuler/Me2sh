#ifndef GEOM_H
#define GEOM_H

#include "gmsh.h"

#include <vector>
#include <array>
#include <iostream>

class Me2sh_Geometry {
    public:
        Me2sh_Geometry(){};
        ~Me2sh_Geometry(){};
        

        std::vector<std::array<double, 3>> points;
        std::vector<std::array<double, 2>> plot_points;
        std::vector<std::array<int, 2>> segments;
        std::vector<int> curveTags;
        std::vector<int> planeTags;
        int firstPointIndex = 0;
        int ExteriorTag;

        void clear(){
            points.clear();
            plot_points.clear();
            segments.clear();
            gmsh::clear();
            firstPointIndex = 0;
        }

        void cleartemp(){
            points.clear();
            plot_points.clear();
            segments.clear();
            firstPointIndex = 0;
        }

        void addSpline(int istart, int iend);
        void addBSpline(int istart, int iend);
#ifdef USE_GEO
        void addBezier(int istart, int iend);
#endif
        void addEllipse(double x, double y, double rx, double ry);
        void addRectangle(double x, double y, double lx, double ly);

        void FuseOverlapping();
        void MakeRectangleAndCut(double x, double y, double lx, double ly);
};

#endif // GEOM_H