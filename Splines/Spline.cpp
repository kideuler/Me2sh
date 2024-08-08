#include "Spline.hpp"

// consider changing to sparse matrices if nv is large
#include <Eigen/Dense> // Eigen dense matrices as computations are relatively small
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

// global 1d quadrature rule
double q[5] = {-(1/3)*sqrt(5+2*sqrt(10/7)), -(1/3)*sqrt(5-2*sqrt(10/7)),0,(1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))};
double w[5] = {(322-13*sqrt(70))/900.0, (322+13*sqrt(70))/900.0, 125/225, (322+13*sqrt(70))/900.0,(322-13*sqrt(70))/900.0};

void Spline2D::init(const Matrix<double> &coords, int degree){
    this->nv = coords.rows;
    this->degree = degree;
    this->xweights = Matrix<double>(nv, degree+1);
    this->yweights = Matrix<double>(nv, degree+1);
    this->params = std::vector<double>(nv+1);
    // copy control points
    control_points.cols = 2;
    for (int i = 0; i<coords.rows; i++){
        control_points.arr.push_back(coords(i,0));
        control_points.arr.push_back(coords(i,1));
        control_points.rows++;
    }

    if (degree == 3){
        Cubic_spline();
    } else {
        std::cerr << "Only degree 3 splines are supported" << std::endl;
    }

    // Gauss quadrature points and weights for arclength computation
    arclength=0.0;

    // reparamaterizing for arc length
    double a,b,I;
    double *temp = new double[nv+1];
    std::array<double,2> xy;
    temp[0] = 0.0;
    for (int i = 0; i<nv; ++i){
        a = params[i];
        b = params[i+1];
        I = 0.0;

        // quadrature
        for (int j = 0; j<5; ++j){
            xy = eval(q[j]*(b-a)/2.0 + (a+b)/2.0);
            I += w[j]*sqrt(xy[0]*xy[0] + xy[1]*xy[1]);
        }
        I = I*(b-a)/2.0;
        arclength += I;
        temp[i+1] = arclength;
    }

    // reparamaterizing
    for (int i = 0; i<nv; ++i){
        params[i] = temp[i]/arclength;
    }
    params[nv] = 1.0;
    delete[] temp;

    int npoints = 1000;
    double dt = 1.0 / (double(npoints));
    double t,x,y;
    double al = 0.0;
    xy = eval(0.0);
    x = xy[0]; y = xy[1];
    for (int i = 1; i<npoints; i++){
        t = i*dt;
        xy = eval(t);
        al += sqrt((xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y));
        x = xy[0]; y = xy[1];
    }
    arclength = al;
}


void Spline2D::Cubic_spline(){
    for (int i = 0; i<nv; ++i){
        params[i] = (double)i / ((double)nv);
    }
    params[nv] = 1.0;

    // setting up Eigen matrices and vectors
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    SpMat A(nv,nv);
    A.reserve(Eigen::VectorXi::Constant(nv,3));
    Eigen::VectorXd bx(nv), by(nv), Dxs(nv), Dys(nv);

    // setting up the tridiagonal matrix and vectors
    for (int ii = 0; ii<nv; ++ii){
        A.insert(ii,ii) = 4.0;
        A.insert(ii,(ii+nv-1)%nv) = 1.0;
        A.insert(ii,(ii+1)%nv) = 4.0;
        bx(ii) = 3.0*(control_points((ii+1)%nv,0) - control_points((ii+nv-1)%nv,0));
        by(ii) = 3.0*(control_points((ii+1)%nv,1) - control_points((ii+nv-1)%nv,1));
    }
    A.makeCompressed();

    // solving the linear system
    solver.analyzePattern(A);
    solver.factorize(A); 
    Dxs = solver.solve(bx); 
    Dys = solver.solve(by); 

    // using solution to set up weights
    for (int ii = 0; ii<nv; ++ii){
        xweights(ii,0) = control_points(ii,0);
        yweights(ii,0) = control_points(ii,1);
        xweights(ii,1) = Dxs(ii);
        yweights(ii,1) = Dys(ii);
        xweights(ii,2) = 3.0*(control_points((ii+1)%nv,0) - control_points(ii,0)) - 2.0*Dxs(ii) - Dxs((ii+1)%nv);
        yweights(ii,2) = 3.0*(control_points((ii+1)%nv,1) - control_points(ii,1)) - 2.0*Dys(ii) - Dys((ii+1)%nv);
        xweights(ii,3) = 2.0*(control_points(ii,0) - control_points((ii+1)%nv,0)) + Dxs(ii) + Dxs((ii+1)%nv);
        yweights(ii,3) = 2.0*(control_points(ii,1) - control_points((ii+1)%nv,1)) + Dys(ii) + Dys((ii+1)%nv);
    }
}

std::array<double,2> Spline2D::eval(double t){

    // checking that t is within 0 to 1
    if (t<0 || t>1){
        std::cerr << "spline evaluation inputs cannot be less than 0 or greater than one: t = " << t << std::endl;
    }

    int n = (int) floor(double(nv)*t);
    if (n==nv){n--;}

    // binary search to find correct segment along spline
    bool stop = false;
    while (!stop) {
        if (n==0){
            // if next param is >= t then stop
            if (params[n+1]>=t){
                stop = true;
            } else {
                n++;
            }
        } else if (n == nv-1){
            if (params[n] <= t){
                stop = true;
            } else {
                n--;
            }
        } else {
            if (params[n] <= t && params[n+1] >= t){
                stop = true;
            } else if (params[n] > t){
                n--;
            } else {
                n++;
            }
        }
    }

    // computing function values
    double x_j = (t-params[n])/(params[n+1]-params[n]);
    double coeff;
    std::array<double,2> xy = {0.0,0.0};
    for (int i = 0; i<degree+1; ++i){
        coeff = 1.0;
        xy[0] += coeff*pow(x_j, i)*xweights(n,i);
        xy[1] += coeff*pow(x_j, i)*yweights(n,i);
    }

    return xy;
}

std::array<double,2> Spline2D::normal(double t){

    // checking that t is within 0 to 1
    if (t<0 || t>1){
        std::cerr << "spline evaluation inputs cannot be less than 0 or greater than one: t = " << t << std::endl;
    }

    int n = (int) floor(double(nv)*t);
    if (n==nv){n--;}

    // binary search to find correct segment along spline
    bool stop = false;
    while (!stop) {
        if (n==0){
            // if next param is >= t then stop
            if (params[n+1]>=t){
                stop = true;
            } else {
                n++;
            }
        } else if (n == nv-1){
            if (params[n] <= t){
                stop = true;
            } else {
                n--;
            }
        } else {
            if (params[n] <= t && params[n+1] >= t){
                stop = true;
            } else if (params[n] > t){
                n--;
            } else {
                n++;
            }
        }
    }

    // computing function values
    double x_j = (t-params[n])/(params[n+1]-params[n]);
    double coeff;
    std::array<double,2> xy = {0.0,0.0};
    for (int i = 1; i<degree+1; ++i){
        coeff = 1.0;
        coeff = coeff*double(i);
        xy[0] += coeff*pow(x_j, i-1)*xweights(n,i);
        xy[1] += coeff*pow(x_j, i-1)*yweights(n,i);
    }

    double nrm = sqrt(xy[0]*xy[0] + xy[1]*xy[1]);
    xy[0] = xy[0]/nrm;
    xy[1] = xy[1]/nrm;
    double temp = xy[0];
    xy[0] = -xy[1];
    xy[1] = temp;

    return xy;
}

void Spline2D::create_segments(double h_target, Matrix<double> &coords, Matrix<int> &segments){
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