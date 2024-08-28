#include "Fem.hpp"

double cs[3][2] = {{2.0/3.0,1.0/6.0},{1.0/6.0,2.0/3.0},{1.0/6.0,1.0/6.0}};
double ws[3] = {1.0/6.0,1.0/6.0,1.0/6.0};

void FemPoisson::assemble(){

    double Jt[4], Jtinv[4];

    int nq = 3;
    int nnodes = 3;
    double tri3[3];
    double tri3_deriv[6] = {-1,-1, 1,0, 0,1};

    double phi[9], dphi[18];
    double grad_sol[6];
    for (int n = 0; n < nq; n++){
        tri3[0] = 1.0 - cs[n][0] - cs[n][1];
        tri3[1] = cs[n][0];
        tri3[2] = cs[n][1];

        phi[nq*n] = tri3[0];
        phi[nq*n+1] = tri3[1];
        phi[nq*n+2] = tri3[2];

        dphi[2*nq*n] = tri3_deriv[0];
        dphi[2*nq*n+1] = tri3_deriv[1];
        dphi[2*nq*n+2] = tri3_deriv[2];
        dphi[2*nq*n+3] = tri3_deriv[3];
        dphi[2*nq*n+4] = tri3_deriv[4];
        dphi[2*nq*n+5] = tri3_deriv[5];
    }

    double ps[6], fs[3];
    double v;
    for (int ii = 0; ii<msh.elems.nrows(); ii++){
            
        for (int nn = 0; nn<nnodes; nn++){
            ps[2*nn] = msh.coords(msh.elems(ii,nn),0);
            ps[2*nn+1] = msh.coords(msh.elems(ii,nn),1);
            fs[nn] = frhs[msh.elems(ii,nn)];
        }
    
        for (int q = 0; q<nq; q++){

            // compute Jacobian
            Jt[0] = 0.0; Jt[1] = 0.0; Jt[2] = 0.0; Jt[3] = 0.0;
            for (int n = 0; n<nnodes; n++){
                for (int jj = 0; jj<2; jj++){
                    for (int kk = 0; kk<2; kk++){
                        Jt[2*jj+kk] += ps[2*n+kk]*dphi[2*nq*q + 2*n + jj];
                    }
                }
            }
            double detJ = (Jt[0]*Jt[3] - Jt[1]*Jt[2]);

            Jtinv[0] = Jt[3]/detJ;
            Jtinv[1] = -Jt[1]/detJ;
            Jtinv[2] = -Jt[2]/detJ;
            Jtinv[3] = Jt[0]/detJ;
            for (int n = 0; n<nnodes; n++){
                grad_sol[2*n] = 0.0;
                grad_sol[2*n+1] = 0.0;
                for (int jj = 0; jj<2; jj++){
                    for (int kk = 0; kk<2; kk++){
                        grad_sol[2*n+jj] += Jtinv[2*jj+kk]*dphi[2*q*nq + 2*n + kk];
                    }
                }
            }

            for (int jj = 0; jj<nnodes; jj++){
                if (!msh.nodes_on_boundary[msh.elems(ii,jj)]){
                    // inserting into right hand side vector
                    double v = 0.0;
                    for (int n = 0; n<nnodes; n++){
                        v += phi[nq*q+n]*fs[n];
                    }
                    b(msh.elems(ii,jj)) += v*ws[q]*detJ*phi[nq*q+jj];

                    // inserting into stiffness matrix
                    for (int kk = 0; kk<nnodes; kk++){
                        for (int dim = 0; dim<2; dim++){
                            A.coeffRef(msh.elems(ii,jj), msh.elems(ii,kk)) += kappa*ws[q]*detJ*grad_sol[2*kk+dim]*grad_sol[2*jj+dim];
                        }
                    }
                }
            }
        }
    }
}

void FemPoisson::apply_dbc(){
    int nbnd = msh.bndnodes.size();
    for (int n = 0; n<msh.coords.nrows(); n++){
        if (msh.nodes_on_boundary[n]){
            A.coeffRef(n,n) = 1.0;
            b(n) = dvals[n];
        }
    }
}

void FemPoisson::solve(){
    // setting up matrix and rhs vector

    A.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(A);
    solver.factorize(A); 
    if (solver.info() == Eigen::Success){
        sol = solver.solve(b); 
    } else {
        std::cout << "Eigen failed to solve" << std::endl;
    }

    for (int n = 0; n<msh.coords.nrows(); n++){
        u[n] = sol(n);
    }
}

void FemPoisson::reset(){
    int nv = 0;
    this->A = Eigen::SparseMatrix<double>(0,0);
    A.reserve(Eigen::VectorXi::Constant(nv,10));
    this->b = Eigen::VectorXd::Zero(nv);
    this->sol = Eigen::VectorXd::Zero(nv);
    this->u = std::vector<double>(nv,0);
}


void FemEikonal::assemble(){

    double Jt[4], Jtinv[4];

    int nq = 3;
    int nnodes = 3;
    double tri3[3];
    double tri3_deriv[6] = {-1,-1, 1,0, 0,1};

    double phi[9], dphi[18];
    double grad_sol[6];
    for (int n = 0; n < nq; n++){
        tri3[0] = 1.0 - cs[n][0] - cs[n][1];
        tri3[1] = cs[n][0];
        tri3[2] = cs[n][1];

        phi[nq*n] = tri3[0];
        phi[nq*n+1] = tri3[1];
        phi[nq*n+2] = tri3[2];

        dphi[2*nq*n] = tri3_deriv[0];
        dphi[2*nq*n+1] = tri3_deriv[1];
        dphi[2*nq*n+2] = tri3_deriv[2];
        dphi[2*nq*n+3] = tri3_deriv[3];
        dphi[2*nq*n+4] = tri3_deriv[4];
        dphi[2*nq*n+5] = tri3_deriv[5];
    }

    double ps[6], fs[3];
    double v;
    for (int ii = 0; ii<msh.elems.nrows(); ii++){
            
        for (int nn = 0; nn<nnodes; nn++){
            ps[2*nn] = msh.coords(msh.elems(ii,nn),0);
            ps[2*nn+1] = msh.coords(msh.elems(ii,nn),1);
        }
    
        for (int q = 0; q<nq; q++){

            // compute Jacobian
            Jt[0] = 0.0; Jt[1] = 0.0; Jt[2] = 0.0; Jt[3] = 0.0;
            for (int n = 0; n<nnodes; n++){
                for (int jj = 0; jj<2; jj++){
                    for (int kk = 0; kk<2; kk++){
                        Jt[2*jj+kk] += ps[2*n+kk]*dphi[2*nq*q + 2*n + jj];
                    }
                }
            }
            double detJ = (Jt[0]*Jt[3] - Jt[1]*Jt[2]);

            Jtinv[0] = Jt[3]/detJ;
            Jtinv[1] = -Jt[1]/detJ;
            Jtinv[2] = -Jt[2]/detJ;
            Jtinv[3] = Jt[0]/detJ;
            for (int n = 0; n<nnodes; n++){
                grad_sol[2*n] = 0.0;
                grad_sol[2*n+1] = 0.0;
                for (int jj = 0; jj<2; jj++){
                    for (int kk = 0; kk<2; kk++){
                        grad_sol[2*n+jj] += Jtinv[2*jj+kk]*dphi[2*q*nq + 2*n + kk];
                    }
                }
            }

            for (int jj = 0; jj<nnodes; jj++){
                if (!msh.nodes_on_boundary[msh.elems(ii,jj)]){

                    double v = 0.0;
                    for (int n = 0; n<nnodes; n++){
                        v += phi[nq*q+n]*1.0;
                    }
                    b(msh.elems(ii,jj)) += v*ws[q]*detJ*phi[nq*q+jj];

                    // inserting into stiffness matrix
                    for (int kk = 0; kk<nnodes; kk++){
                        for (int dim = 0; dim<2; dim++){
                            A.coeffRef(msh.elems(ii,jj), msh.elems(ii,kk)) += 1*ws[q]*detJ*grad_sol[2*kk+dim]*grad_sol[2*jj+dim];
                        }
                        //A.coeffRef(msh.elems(ii,jj), msh.elems(ii,kk)) += -alpha*alpha*ws[q]*detJ*phi[nq*q+jj]*phi[nq*q+kk];
                    }
                }
            }
        }
    }
}

void FemEikonal::apply_dbc(){
    int nbnd = msh.bndnodes.size();
    for (int n = 0; n<msh.coords.nrows(); n++){
        if (msh.nodes_on_boundary[n]){
            A.coeffRef(n,n) = 1.0;
            b(n) = 0.0;
        }
    }
}

void FemEikonal::solve(){
    // setting up matrix and rhs vector

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

    for (int n = 0; n<msh.coords.nrows(); n++){
        u[n] = sol(n);
    }


}


void FemEikonal::reset(){
    int nv = 0;
    this->A = Eigen::SparseMatrix<double>(0,0);
    A.reserve(Eigen::VectorXi::Constant(nv,10));
    this->b = Eigen::VectorXd::Zero(nv);
    this->sol = Eigen::VectorXd::Zero(nv);
    this->u = std::vector<double>(nv,0);
}