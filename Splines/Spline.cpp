#include "Spline.hpp"

// consider changing to sparse matrices if nv is large
#include <Eigen/Dense> // Eigen dense matrices as computations are relatively small

void Spline2D::init(Matrix<double> coords, int degree){
    this->nv = coords.rows;
    this->degree = degree;
    this->control_points = coords;
    this->xweights = Matrix<double>(nv, degree+1);
    this->yweights = Matrix<double>(nv, degree+1);
    this->params = std::vector<double>(nv+1);

    if (degree == 3){
        Cubic_spline();
    } else {
        std::cerr << "Only degree 3 splines are supported" << std::endl;
    }

    // Gauss quadrature points and weights for arclength computation
    double arclength=0.0;
    double q[5] = {-(1/3)*sqrt(5+2*sqrt(10/7)), -(1/3)*sqrt(5-2*sqrt(10/7)),0,(1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))};
    double w[5] = {(322-13*sqrt(70))/900.0, (322+13*sqrt(70))/900.0, 125/225, (322+13*sqrt(70))/900.0,(322-13*sqrt(70))/900.0};

    // reparamaterizing for arc length
    
}


void Spline2D::Cubic_spline(){

    for (int i = 0; i<nv; ++i){
        params[i] = (double)i / ((double)nv+1);
    }
    params[nv] = 1.0;

    // setting up Eigen matrices and vectors
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nv, nv);
    Eigen::VectorXd bx(nv), by(nv), Dxs(nv), Dys(nv);

    // setting up the tridiagonal matrix and vectors
    for (int ii = 0; ii<nv; ++ii){
        A(ii,ii) = 4.0;
        A(ii,(ii+nv-1)%nv) = 1.0;
        A(ii,(ii+1)%nv) = 4.0;
        bx(ii) = 3.0*(control_points((ii+1)%nv,0) - control_points((ii+nv-1)%nv,0));
        by(ii) = 3.0*(control_points((ii+1)%nv,1) - control_points((ii+nv-1)%nv,1));
    }

    // solving the linear system
    Dxs = A.lu().solve(bx);
    Dys = A.lu().solve(by);

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