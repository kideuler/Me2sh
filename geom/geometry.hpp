#ifndef GEOM_H
#define GEOM_H

#include "gmsh.h"

#include <vector>
#include <array>
#include <iostream>

class Me2sh_Geometry {
    public:
        Me2sh_Geometry(){if (gmsh::isInitialized() == 0){gmsh::initialize(); gmsh::model::add("me2sh"); inititalized_in_geom = true;}};
        ~Me2sh_Geometry(){if (inititalized_in_geom){gmsh::finalize();}};
        

        std::vector<std::array<double, 3>> points;
        std::vector<std::array<double, 2>> plot_points;
        std::vector<std::array<int, 2>> segments;
        std::vector<int> curveTags;
        std::vector<int> planeTags;
        int firstPointIndex = 0;
        int ExteriorTag = -1;

        void clear(){
            points.clear();
            plot_points.clear();
            segments.clear();
            curveTags.clear();
            planeTags.clear();
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

        void FuseOverlapping(bool exterior = false);
        void MakeRectangleAndCut(double x, double y, double lx, double ly);


    private:
        bool inititalized_in_geom = false;
};

#endif // GEOM_H