#include "Bezier.hpp"

void Bezier2D::init(Matrix<double> coords){
    control_points = coords;
    control_points.cols = 2;
    control_points.arr.push_back(coords(0,0));
    control_points.arr.push_back(coords(0,1));
    control_points.rows++;
    nv = control_points.rows;
    arclength = 0.0;

    int npoints = 1000;
    double dt = 1.0 / (double(npoints));
    double t,x,y;
    double al = 0.0;
    std::array<double,2> xy = eval(0.0);
    x = xy[0]; y = xy[1];
    for (int i = 1; i<npoints; i++){
        t = i*dt;
        xy = eval(t);
        al += sqrt((xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y));
        x = xy[0]; y = xy[1];
    }
    arclength = al;
}

std::array<double,2> Bezier2D::eval(double t){
    std::array<double,2> xy;
    
    Matrix<double> B = control_points;
    for (int i = 1; i<=nv; i++){
        for (int j = 0; j<(nv-i); j++){
            B(j,0) = (1-t)*B(j,0) + t*B(j+1,0);
            B(j,1) = (1-t)*B(j,1) + t*B(j+1,1);
        }
    }
    xy[0] = B(0,0);
    xy[1] = B(0,1);
    return xy;
}

std::array<double,2> Bezier2D::normal(double t){
    std::array<double,2> xy;
    Matrix<double> B = control_points;
    Matrix<double> C = control_points;
    for (int i = 1; i<=nv-1; i++){
        for (int j = 0; j<(nv-i); j++){
            B(j,0) = (1-t)*B(j,0) + t*B(j+1,0);
            B(j,1) = (1-t)*B(j,1) + t*B(j+1,1);
        }
    }

    for (int i = 2; i<=nv; i++){
        for (int j = 1; j<=(nv-i); j++){
            C(j,0) = (1-t)*C(j,0) + t*C(j+1,0);
            C(j,1) = (1-t)*C(j,1) + t*C(j+1,1);
        }
    }
    xy[0] = double(nv)*(C(1,0) - B(0,0));
    xy[1] = double(nv)*(C(1,1) - B(0,1));
    double nrm = sqrt(xy[0]*xy[0] + xy[1]*xy[1]);
    xy[0] = xy[0]/nrm;
    xy[1] = xy[1]/nrm;
    double temp = xy[0];
    xy[0] = -xy[1];
    xy[1] = temp;
    return xy;
}

void Bezier2D::create_segments(double h_target, Matrix<double> &coords, Matrix<int> &segments){
    int npoints = (int) (arclength/(h_target));
    double dt = 1.0 / (double(npoints));

    coords = Matrix<double>(0,0);
    coords.cols = 2;
    coords.rows = npoints;
    coords.arr.resize(npoints*2);

    segments = Matrix<int>(0,0);
    segments.cols = 2;
    segments.rows = npoints;
    segments.arr.resize(npoints*2);

    std::array<double,2> xy;
    double t;
    for (double i = 0; i<npoints; i++){
        t = i*dt;
        xy = eval(t);
        coords(i,0) = xy[0];
        coords(i,1) = xy[1];

        segments(i,0) = i;
        segments(i,1) = ((int)i+1)%npoints;
    }
}