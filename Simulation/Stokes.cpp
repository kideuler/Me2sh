#include "Stokes.hpp"
#include<Eigen/IterativeLinearSolvers>

double cs[3][2] = {{2.0/3.0,1.0/6.0},{1.0/6.0,2.0/3.0},{1.0/6.0,1.0/6.0}};
double ws[3] = {1.0/6.0,1.0/6.0,1.0/6.0};

static inline void Bubble_shapefunction(double x, double y, double *phi, double *dphi){
    phi[0] = 9*x*x*y + 9*x*y*y - 9*x*y - x - y + 1;
    phi[1] = x*(9*x*y + 9*y*y - 9*y + 1);
    phi[2] = y*(9*x*x + 9*x*y - 9*x + 1);
    phi[3] = 27*x*y*(-x-y+1);


    dphi[0] = 18*x*y + 9*y*y - 9*y - 1;
    dphi[1] = 18*x*y - 9*x + 9*x*x - 1;

    dphi[2] = 18*x*y + 9*y*y - 9*y + 1;
    dphi[3] = x*(9*x + 18*y - 9);

    dphi[4] = y*(18*x + 9*y - 9);
    dphi[5] = 9*x*x + 18*x*y - 9*x + 1;

    dphi[6] = -27*y*(2*x + y - 1);
    dphi[7] = -27*x*(x + 2*y - 1);
}

static inline void Langrange_shapefunction(double x, double y, double *phi, double *dphi){
    phi[0] = 1 - x - y;
    phi[1] = x;
    phi[2] = y;

    dphi[0] = -1;
    dphi[1] = -1;
    dphi[2] = 1;
    dphi[3] = 0;
    dphi[4] = 0;
    dphi[5] = 1;
}


// assembly using Bubble function for velocity and Lagrange for pressure
void FemStokes::assemble_B1_P1(){

    // Jacobians and inverses
    double Jt_P1[4], Jtinv_P1[4];
    double Jt_B1[4], Jtinv_B1[4];

    // number of quadrature points
    int nq = 3;
    int nnodes = 3;
    int offset = nv + nelems;
    std::cout << "offset: " << offset << std::endl;

    // shape function and its derivatives
    double P1_phi[nq][3], P1_dphi[nq][6];
    double P1_grad_sol[6];
    double B1_phi[nq][4], B1_dphi[nq][8];
    double B1_grad_sol[8];
    double B1_div_sol[4];

    // compute shape functions
    for (int n = 0; n < nq; n++){
        double x = cs[n][0];
        double y = cs[n][1];

        // Lagrange shape functions
        Langrange_shapefunction(x, y, P1_phi[n], P1_dphi[n]);

        // Bubble shape functions
        Bubble_shapefunction(x, y, B1_phi[n], B1_dphi[n]);
    }

    double ps[8]; // coordinates of triangle and center
    double fs_x[4], fs_y[4]; // right hand side vector
    double v_x,v_y;
    int elemnodes[4]; // 3 nodes + 1 center node
    
    // loop through all elements and construct elemental stiffness matrix
    for (int ii = 0; ii<nelems; ii++){
        
        // get coordinates of triangle and add center node, also get rhs values
        for (int nn = 0; nn<nnodes; nn++){
            elemnodes[nn] = msh.elems(ii,nn);
            ps[2*nn] = msh.coords(msh.elems(ii,nn),0);
            ps[2*nn+1] = msh.coords(msh.elems(ii,nn),1);
            fs_x[nn] = frhs_x[msh.elems(ii,nn)];
            fs_y[nn] = frhs_y[msh.elems(ii,nn)];
        }
        elemnodes[3] = nv + ii; // center node
        ps[6] = msh.coords(elemnodes[3],0);
        ps[7] = msh.coords(elemnodes[3],1);
        fs_x[3] = frhs_x[elemnodes[3]];
        fs_y[3] = frhs_y[elemnodes[3]];

        // loop over quadrature points
        for (int q = 0; q<nq; q++){

            // compute Jacobian P1
            Jt_P1[0] = 0.0; Jt_P1[1] = 0.0; Jt_P1[2] = 0.0; Jt_P1[3] = 0.0;
            for (int n = 0; n<3; n++){
                for (int jj = 0; jj<2; jj++){
                    for (int kk = 0; kk<2; kk++){
                        Jt_P1[2*jj+kk] += ps[2*n+kk]*P1_dphi[q][2*n + jj];
                    }
                }
            }
            double detJ_P1 = (Jt_P1[0]*Jt_P1[3] - Jt_P1[1]*Jt_P1[2]);
            Jtinv_P1[0] = Jt_P1[3]/detJ_P1;
            Jtinv_P1[1] = -Jt_P1[1]/detJ_P1;
            Jtinv_P1[2] = -Jt_P1[2]/detJ_P1;
            Jtinv_P1[3] = Jt_P1[0]/detJ_P1;


            // compute Jacobian B1
            Jt_B1[0] = 0.0; Jt_B1[1] = 0.0; Jt_B1[2] = 0.0; Jt_B1[3] = 0.0;
            for (int n = 0; n<4; n++){
                for (int jj = 0; jj<2; jj++){
                    for (int kk = 0; kk<2; kk++){
                        Jt_B1[2*jj+kk] += ps[2*n+kk]*B1_dphi[q][2*n + jj];
                    }
                }
            }
            double detJ_B1 = (Jt_B1[0]*Jt_B1[3] - Jt_B1[1]*Jt_B1[2]);
            Jtinv_B1[0] = Jt_B1[3]/detJ_B1;
            Jtinv_B1[1] = -Jt_B1[1]/detJ_B1;
            Jtinv_B1[2] = -Jt_B1[2]/detJ_B1;
            Jtinv_B1[3] = Jt_B1[0]/detJ_B1;

            if (detJ_B1 < 0.0 || detJ_P1 < 0.0){
                std::cout << "detJ_B1: " << detJ_B1 << " detJ_P1: " << detJ_P1 << std::endl;
                continue;
            }


            // compute P1 gradient
            for (int n = 0; n<3; n++){
                P1_grad_sol[2*n] = 0.0;
                P1_grad_sol[2*n+1] = 0.0;
                for (int jj = 0; jj<2; jj++){
                    for (int kk = 0; kk<2; kk++){
                        P1_grad_sol[2*n+jj] += Jtinv_P1[2*jj+kk]*P1_dphi[q][2*n+kk];
                    }
                }
            }

            // compute B1 gradient and divergence
            for (int n = 0; n<4; n++){
                B1_grad_sol[2*n] = 0.0;
                B1_grad_sol[2*n+1] = 0.0;
                for (int jj = 0; jj<2; jj++){
                    for (int kk = 0; kk<2; kk++){
                        B1_grad_sol[2*n+jj] += Jtinv_B1[2*jj+kk]*B1_dphi[q][2*n+kk];
                    }
                }
                B1_div_sol[n] = B1_grad_sol[2*n] + B1_grad_sol[2*n+1];
            }

            // inserting B1 stiffness matrices and rhs vectors
            for (int jj = 0; jj<4; jj++){
                //if (jj != 3 && msh.nodes_on_boundary[elemnodes[jj]]){continue;}

                // inserting mu*grad(u)*grad(v) into matrix
                for (int kk = 0; kk<4; kk++){
                    for (int dim = 0; dim<2; dim++){
                        A.coeffRef(elemnodes[jj],elemnodes[kk]) += mu*ws[q]*detJ_B1*B1_grad_sol[2*kk+dim]*B1_grad_sol[2*jj+dim];
                        A.coeffRef(elemnodes[jj]+offset,elemnodes[kk]+offset) += mu*ws[q]*detJ_B1*B1_grad_sol[2*kk+dim]*B1_grad_sol[2*jj+dim];
                    }
                }

                // inserting -p*div(v) into matrix
                for (int kk = 0; kk<3; kk++){
                    // top right side
                    A.coeffRef(elemnodes[jj],elemnodes[kk]+2*offset) += -ws[q]*detJ_B1*P1_phi[q][kk]*B1_div_sol[jj];
                    A.coeffRef(elemnodes[jj]+offset,elemnodes[kk]+2*offset) += -ws[q]*detJ_B1*P1_phi[q][kk]*B1_div_sol[jj];
                }


                // inserting into right hand side vector
                v_x = 0.0; v_y = 0.0;
                for (int kk = 0; kk<4; kk++){
                    v_x += B1_phi[q][kk]*fs_x[kk];
                    v_y += B1_phi[q][kk]*fs_y[kk];
                }
                b(elemnodes[jj]) += v_x*ws[q]*detJ_B1*B1_phi[q][jj];
                b(elemnodes[jj]+offset) += v_y*ws[q]*detJ_B1*B1_phi[q][jj];
            }


            // inserting pressure terms
            for (int jj = 0; jj<3; jj++){
                //if (msh.nodes_on_boundary[elemnodes[jj]]){continue;}
                
                // inserting -eps*p*q into matrix
                for (int kk = 0; kk<3; kk++){
                    A.coeffRef(elemnodes[jj]+2*offset,elemnodes[kk]+2*offset) += eps*ws[q]*detJ_P1*P1_phi[q][jj]*P1_phi[q][kk];
                }
                b(elemnodes[jj]+2*offset) = eps;

                // inserting -q*div(u) into matrix
                for (int kk = 0; kk<4; kk++){
                    A.coeffRef(elemnodes[jj]+2*offset,elemnodes[kk]) += ws[q]*detJ_P1*B1_div_sol[kk]*P1_phi[q][jj];
                    A.coeffRef(elemnodes[jj]+2*offset,elemnodes[kk]+offset) += ws[q]*detJ_P1*B1_div_sol[kk]*P1_phi[q][jj];
                }

            }
        }
    }
}

int E = 1e6;
void FemStokes::apply_dbc(){
    // apply Dirichlet boundary conditions
    int offset = nv + nelems;
    for (int ii = 0; ii < nv; ii++){
        if (msh.nodes_on_boundary[ii]){
            // Zero out the rest of the row for x-component
            // for (Eigen::SparseMatrix<double>::InnerIterator it(A, ii); it; ++it) {
            //     if (it.row() != ii) {
            //         A.coeffRef(ii, it.col()) = 0.0;
            //     }
            // }
            // Set velocity to zero for x-component
            A.coeffRef(ii, ii) += E;
            b(ii) += E*dvals_x[ii];

            // Zero out the rest of the row for y-component
            // for (Eigen::SparseMatrix<double>::InnerIterator it(A, ii + offset); it; ++it) {
            //     if (it.row() != ii + offset) {
            //         A.coeffRef(ii + offset, it.col()) = 0.0;
            //     }
            // }
            // Set velocity to zero for y-component
            A.coeffRef(ii + offset, ii + offset) += E;
            b(ii + offset) += E*dvals_y[ii];
        }
    }
}

void FemStokes::solve(){
    // setting up matrix and rhs vector
    int offset = nv + nelems;
    A.makeCompressed();

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    sol = solver.solve(b); 
    if (solver.info() == Eigen::Success){
        
        std::cout << "#iterations:     " << solver.iterations() << std::endl;
        std::cout << "estimated error: " << solver.error()      << std::endl;
    } else {
        std::cout << "Eigen failed to solve" << std::endl;
    }

    // A.makeCompressed();
    // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    // solver.analyzePattern(A);
    // solver.factorize(A); 
    // if (solver.info() == Eigen::Success){
    //     sol = solver.solve(b); 
    // } else {
    //     std::cout << "Eigen failed to solve" << std::endl;
    // }

    // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > lscg;
    // lscg.compute(A);
    // sol = lscg.solve(b);
    // std::cout << "#iterations:     " << lscg.iterations() << std::endl;
    // std::cout << "estimated error: " << lscg.error()      << std::endl;
    // // update b, and solve again
    // sol = lscg.solve(b);

    for (int ii = 0; ii<nv; ii++){
        //std::cout << "sol: " << sol(ii+2*offset) << std::endl;
        u_x[ii] = sol(ii);
        u_y[ii] = sol(ii + offset);
        p[ii] = sol(ii + 2*offset);
    }
}