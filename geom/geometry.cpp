#include "geometry.hpp"

void Me2sh_Geometry::addEllipse(double x, double y, double rx, double ry){
    int tag = gmsh::model::occ::addEllipse(x, y, 0.0, rx, ry);
    gmsh::model::occ::synchronize();
    curveTags.push_back(tag);

    int wtag = gmsh::model::occ::addCurveLoop({tag});
    gmsh::model::occ::synchronize();

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

void Me2sh_Geometry::addRectangle(double x, double y, double lx, double ly){
    int v1 = gmsh::model::occ::addPoint(x-lx, y-ly, 0.0);
    int v2 = gmsh::model::occ::addPoint(x+lx, y-ly, 0.0);
    int v3 = gmsh::model::occ::addPoint(x+lx, y+ly, 0.0);
    int v4 = gmsh::model::occ::addPoint(x-lx, y+ly, 0.0);

    std::vector<int> ltags(4);
    ltags[0] = gmsh::model::occ::addLine(v1, v2);
    curveTags.push_back(ltags[0]);
    ltags[1] = gmsh::model::occ::addLine(v2, v3);
    curveTags.push_back(ltags[1]);
    ltags[2] = gmsh::model::occ::addLine(v3, v4);
    curveTags.push_back(ltags[2]);
    ltags[3] = gmsh::model::occ::addLine(v4, v1);
    curveTags.push_back(ltags[3]);

    gmsh::model::occ::synchronize();

    int wtag = gmsh::model::occ::addCurveLoop(ltags);
    gmsh::model::occ::synchronize();

    int planetag = gmsh::model::occ::addPlaneSurface({wtag});
    gmsh::model::occ::synchronize();

    planeTags.push_back(planetag);

    plot_points.clear();
    for (int i = 0; i < 4; i++){
        std::vector<double> min, max;
        gmsh::model::getParametrizationBounds(1, ltags[i], min, max);

        std::vector<double> params;
        double h = (max[0]-min[0])/20.0;
        for (double t = min[0]; t < max[0]; t += h){
            params.push_back(t);
        }
        std::vector<double> coords;

        gmsh::model::getValue(1, ltags[i], params, coords);
        for (int i = 0; i < coords.size(); i+=3){
            plot_points.push_back({coords[i], coords[i+1]});
        }
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
    curveTags.push_back(stag);

    int wtag = gmsh::model::occ::addCurveLoop({stag});
    gmsh::model::occ::synchronize();

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
    curveTags.push_back(stag);

    int wtag = gmsh::model::occ::addCurveLoop({stag});
    gmsh::model::occ::synchronize();

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

#ifdef USE_GEO
void Me2sh_Geometry::addBezier(int istart, int iend){

    std::vector<int> tags;

    int gstart = gmsh::model::geo::addPoint(points[istart][0], points[istart][1], points[istart][2]);
    tags.push_back(gstart);
    for (int i = istart+1; i < iend-1; i++){
        int tag = gmsh::model::geo::addPoint(points[i][0], points[i][1], points[i][2]);
        tags.push_back(tag);
    }
    int gend = gmsh::model::geo::addPoint(points[iend-1][0], points[iend-1][1], points[iend-1][2]);
    gmsh::model::geo::synchronize();
    tags.push_back(gend);
    tags.push_back(gstart);
    

    firstPointIndex = iend;

    int stag = gmsh::model::geo::addBezier(tags);
    gmsh::model::geo::synchronize();
    curveTags.push_back(stag);

    int wtag = gmsh::model::geo::addCurveLoop({stag});
    gmsh::model::geo::synchronize();

    int planetag = gmsh::model::geo::addPlaneSurface({wtag});
    gmsh::model::geo::synchronize();

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
#endif