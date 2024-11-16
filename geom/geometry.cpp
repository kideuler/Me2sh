#include "geometry.hpp"

void Me2sh_Geometry::addEllipse(double x, double y, double rx, double ry){
    int tag = gmsh::model::occ::addEllipse(x, y, 0.0, rx, ry);
    gmsh::model::occ::synchronize();

    int wtag = gmsh::model::occ::addCurveLoop({tag});
    gmsh::model::occ::synchronize();
    curveLoopTags.push_back(wtag);

    int planetag = gmsh::model::occ::addPlaneSurface({wtag});
    gmsh::model::occ::synchronize();

    planeTags.push_back(planetag);

    std::vector<double> min, max;
    gmsh::model::getParametrizationBounds(1, tag, min, max);

    std::vector<double> params;
    double h = (max[0]-min[0])/100.0;
    for (double t = min[0]; t < max[0]; t += h){
        params.push_back(t);
    }
    std::vector<double> coords;

    gmsh::model::getValue(1,tag,params,coords);
    plot_points.clear();
    for (int i = 0; i < coords.size(); i+=3){
        plot_points.push_back({coords[i], coords[i+1]});
    }
}

void Me2sh_Geometry::addSpline(int istart, int iend){

    std::vector<int> tags;

    int gstart = gmsh::model::occ::addPoint(points[istart][0], points[istart][1], points[istart][2]);
    tags.push_back(gstart);
    for (int i = istart+1; i < iend-1; i++){
        int tag = gmsh::model::occ::addPoint(points[i][0], points[i][1], points[i][2]);
        tags.push_back(tag);
    }
    int gend = gmsh::model::occ::addPoint(points[iend-1][0], points[iend-1][1], points[iend-1][2]);
    tags.push_back(gend);
    tags.push_back(gstart);

    firstPointIndex = iend;

    int stag = gmsh::model::occ::addSpline(tags);
    gmsh::model::occ::synchronize();

    int wtag = gmsh::model::occ::addCurveLoop({stag});
    gmsh::model::occ::synchronize();
    curveLoopTags.push_back(wtag);

    int planetag = gmsh::model::occ::addPlaneSurface({wtag});
    gmsh::model::occ::synchronize();

    planeTags.push_back(planetag);

    std::vector<double> min, max;
    gmsh::model::getParametrizationBounds(1, stag, min, max);

    std::vector<double> params;
    double h = (max[0]-min[0])/500.0;
    for (double t = min[0]; t < max[0]; t += h){
        params.push_back(t);
    }
    std::vector<double> coords;

    gmsh::model::getValue(1,stag,params,coords);
    plot_points.clear();
    for (int i = 0; i < coords.size(); i+=3){
        plot_points.push_back({coords[i], coords[i+1]});
    }
}

void Me2sh_Geometry::addBSpline(int istart, int iend){

    std::vector<int> tags;

    int gstart = gmsh::model::occ::addPoint(points[istart][0], points[istart][1], points[istart][2]);
    tags.push_back(gstart);
    for (int i = istart+1; i < iend-1; i++){
        int tag = gmsh::model::occ::addPoint(points[i][0], points[i][1], points[i][2]);
        tags.push_back(tag);
    }
    int gend = gmsh::model::occ::addPoint(points[iend-1][0], points[iend-1][1], points[iend-1][2]);
    tags.push_back(gend);
    tags.push_back(gstart);

    firstPointIndex = iend;

    int stag = gmsh::model::occ::addBSpline(tags);
    gmsh::model::occ::synchronize();

    int wtag = gmsh::model::occ::addCurveLoop({stag});
    gmsh::model::occ::synchronize();
    curveLoopTags.push_back(wtag);

    int planetag = gmsh::model::occ::addPlaneSurface({wtag});
    gmsh::model::occ::synchronize();

    planeTags.push_back(planetag);

    std::vector<double> min, max;
    gmsh::model::getParametrizationBounds(1, stag, min, max);

    std::vector<double> params;
    double h = (max[0]-min[0])/500.0;
    for (double t = min[0]; t < max[0]; t += h){
        params.push_back(t);
    }
    std::vector<double> coords;

    gmsh::model::getValue(1,stag,params,coords);
    plot_points.clear();
    for (int i = 0; i < coords.size(); i+=3){
        plot_points.push_back({coords[i], coords[i+1]});
    }
}

void Me2sh_Geometry::addBezier(int istart, int iend){

    std::vector<int> tags;

    int gstart = gmsh::model::occ::addPoint(points[istart][0], points[istart][1], points[istart][2]);
    tags.push_back(gstart);
    for (int i = istart+1; i < iend-1; i++){
        int tag = gmsh::model::occ::addPoint(points[i][0], points[i][1], points[i][2]);
        tags.push_back(tag);
    }
    int gend = gmsh::model::occ::addPoint(points[iend-1][0], points[iend-1][1], points[iend-1][2]);
    gmsh::model::occ::synchronize();
    tags.push_back(gend);
    tags.push_back(gstart);
    

    firstPointIndex = iend;

    int stag = gmsh::model::occ::addBezier(tags);
    gmsh::model::occ::synchronize();

    int wtag = gmsh::model::occ::addCurveLoop({stag});
    gmsh::model::occ::synchronize();
    curveLoopTags.push_back(wtag);

    int planetag = gmsh::model::occ::addPlaneSurface({wtag});
    gmsh::model::occ::synchronize();

    planeTags.push_back(planetag);

    std::vector<double> min, max;
    gmsh::model::getParametrizationBounds(1, stag, min, max);

    std::vector<double> params;
    double h = (max[0]-min[0])/500.0;
    for (double t = min[0]; t < max[0]; t += h){
        params.push_back(t);
    }
    std::vector<double> coords;

    gmsh::model::getValue(1,stag,params,coords);
    plot_points.clear();
    for (int i = 0; i < coords.size(); i+=3){
        plot_points.push_back({coords[i], coords[i+1]});
    }
}