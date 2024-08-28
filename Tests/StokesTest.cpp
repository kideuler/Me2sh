#include "Stokes.hpp"

double mu = 1.0;

std::vector<double> vel_x(Matrix<double> xs){
    int nv = xs.nrows();
    std::vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = xs(i,0)*xs(i,0) + xs(i,1)*xs(i,1);
    }
    return u;
}

std::vector<double> vel_y(Matrix<double> xs){
    int nv = xs.nrows();
    std::vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = 2*xs(i,0)*xs(i,0) - 2*xs(i,0)*xs(i,1);
    }
    return u;
}

std::vector<double> pexpr(Matrix<double> xs){
    int nv = xs.nrows();
    std::vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = xs(i,0) + xs(i,1) - 1;
    }
    return u;
}

std::vector<double> frhs_x(Matrix<double> xs){
    int nv = xs.nrows();
    std::vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = 1-4*mu;
    }
    return u;
}

std::vector<double> frhs_y(Matrix<double> xs){
    int nv = xs.nrows();
    std::vector<double> u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = 1-4*mu;
    }
    return u;
}


int main(){

    int npoints = 20;
    IntMatrix segments;
    DoubleMatrix coords;
    double h = Circle(&segments, &coords, npoints, true);
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
    msh.Refine(0.05);
    msh.Smooth(100);
    msh.compute_boundary_nodes();

    // print number of nodes and elements
    std::cout << "Number of nodes: " << msh.coords.nrows() << std::endl;
    std::cout << "Number of elements: " << msh.elems.nrows() << std::endl;

    // initializing the Stokes solver class
    FemStokes fem;
    fem.init(msh, mu, frhs_x(msh.coords), frhs_y(msh.coords), vel_x(msh.coords), vel_y(msh.coords));
    fem.make_3P1_mesh();
    fem.dvals_x = vel_x(fem.msh.coords);
    fem.dvals_y = vel_y(fem.msh.coords);
    fem.frhs_x = frhs_x(fem.msh.coords);
    fem.frhs_y = frhs_y(fem.msh.coords);
    fem.assemble_B1_P1();
    fem.apply_dbc();
    fem.solve();

    int nv = msh.coords.nrows();
    std::vector<double> error_ux(nv), error_uy(nv), error_p(nv);
    double maxE_ux = 0.0;
    double maxE_uy = 0.0;
    double maxE_p = 0.0;
    double l2E_ux = 0.0;
    double l2E_uy = 0.0;
    double l2E_p = 0.0;

    std::vector<double> u_x = vel_x(msh.coords);
    std::vector<double> u_y = vel_y(msh.coords);
    std::vector<double> p = pexpr(msh.coords);

    for (int n = 0; n<nv; n++) {
        error_ux[n] = abs(u_x[n]-fem.u_x[n]); 
        error_uy[n] = abs(u_y[n]-fem.u_y[n]); 
        error_p[n] = abs(p[n]-fem.p[n]); 
        // std::cout << "u_x: " << u_x[n] << " fem.u_x: " << fem.u_x[n] << " error_ux: " << error_ux[n] << std::endl;
        // std::cout << "u_y: " << u_y[n] << " fem.u_y: " << fem.u_y[n] << " error_uy: " << error_uy[n] << std::endl;
        // std::cout << "p: " << p[n] << " fem.p: " << fem.p[n] << " error_p: " << error_p[n] << std::endl;

        l2E_ux += error_ux[n]*error_ux[n];
        l2E_uy += error_uy[n]*error_uy[n];
        l2E_p += error_p[n]*error_p[n];

        if(error_ux[n]>maxE_ux){maxE_ux=error_ux[n];}
        if(error_uy[n]>maxE_uy){maxE_uy=error_uy[n];}
        if(error_p[n]>maxE_p){maxE_p=error_p[n];}
    }
    l2E_ux = sqrt(l2E_ux);
    l2E_uy = sqrt(l2E_uy);
    l2E_p = sqrt(l2E_p);

    std::cout << "l_inf error u_x: " << maxE_ux << std::endl;
    std::cout << "l_inf error u_y: " << maxE_uy << std::endl;
    std::cout << "l_inf error p: " << maxE_p << std::endl;

    std::cout << "l_2 error u_x: " << l2E_ux << std::endl;
    std::cout << "l_2 error u_y: " << l2E_uy << std::endl;
    std::cout << "l_2 error p: " << l2E_p << std::endl;

    if (l2E_ux>1e-2 || maxE_ux>1e-2 || l2E_uy>1e-2 || maxE_uy>1e-2 || l2E_p>1e-2 || maxE_p>1e-2){
        std::cout << "Test failed" << std::endl;
        return 0;
    } else {
        std::cout << "Test passed" << std::endl;
        return 1;
    }

    

    return 0;
}