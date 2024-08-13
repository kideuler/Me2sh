#include "Fem.hpp"

std::vector<double> uexpr(Matrix<double> xs){
    int nv = xs.nrows();
    std::vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = 10*sin(2*xs(i,0))*cos(2*xs(i,1));
    }
    return u;
}

std::vector<double> uexpr2(Matrix<double> xs){
    int nv = xs.nrows();
    std::vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = xs(i,0)*xs(i,0) + xs(i,1)*xs(i,1);
    }
    return u;
}

std::vector<double> ulap(Matrix<double> xs, double kappa){
    int nv = xs.nrows();
    std::vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = -80*sin(2*xs(i,0))*cos(2*xs(i,1))*(-kappa);
    }
    return u;
}

std::vector<double> ulap2(Matrix<double> xs, double kappa){
    int nv = xs.nrows();
    std::vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = -4;
    }
    return u;
}


int TEST1(){

    int npoints = 100;
    IntMatrix segments;
    DoubleMatrix coords;
    double h = Circle(&segments, &coords, npoints, false);
    Matrix<double> coords_msh(npoints,2);
    Matrix<int> segments_msh(npoints,2);

    for (int i = 0; i<npoints; i++){
        coords_msh(i,0) = coords.data[2*i];
        coords_msh(i,1) = coords.data[2*i+1];
        segments_msh(i,0) = segments.data[2*i];
        segments_msh(i,1) = segments.data[2*i+1];
    }

    Mesh2D msh;

    // triangulating cooreinates
    msh.Triangulate(coords_msh,segments_msh);
    msh.Refine(0.01);
    msh.Smooth(100);
    msh.compute_boundary_nodes();

    double kappa = 1.0;
    std::vector<double> u = uexpr(msh.coords); // true solution
    std::vector<double> frhs = ulap(msh.coords, kappa); // rhs expression

    // initializing the Poisson solver class
    FemPoisson fem;
    fem.init(msh, frhs, kappa, u);
    fem.assemble();
    fem.apply_dbc();
    fem.solve();
    
    int nv = msh.coords.nrows();
    std::vector<double> error(nv);
    double maxE = 0.0;
    double l2E = 0.0;
    for (int n = 0; n<nv; n++) {
        error[n] = abs(u[n]-fem.u[n]); 
        l2E += error[n]*error[n];
        if(error[n]>maxE){maxE=error[n];}
    }
    l2E = sqrt(l2E);
    std::cout << "l_inf error: " << maxE << std::endl;
    std::cout << "l_2 error: " << l2E << std::endl;
    if (l2E>1e-2){
        std::cout << "Test failed" << std::endl;
        return 0;
    } else {
        std::cout << "Test passed" << std::endl;
        return 1;
    }
}

int main(){
    if (TEST1()){
        return 0;
    } else {
        return 1;
    }
}