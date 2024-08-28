#ifndef STOKES_HPP
#define STOKES_HPP

#include <Eigen/Sparse>
#include "TriMesh.hpp"
#include <math.h>

class FemStokes {

    public:

        FemStokes(){};

        void make_3P1_mesh(){
            for (int ii = 0; ii<nelems; ii++){
                double xc = 0.0;
                double yc = 0.0;
                for (int nn = 0; nn<3; nn++){
                    xc += this->msh.coords(msh.elems(ii,nn),0);
                    yc += this->msh.coords(msh.elems(ii,nn),1);
                }
                xc /= 3.0;
                yc /= 3.0;

                this->msh.coords.arr.push_back(xc);
                this->msh.coords.arr.push_back(yc);
                this->msh.coords.rows++;
            }
        }

        void init(Mesh2D &msh, double mu, const std::vector<double> &frhs_x, const std::vector<double> &frhs_y, const std::vector<double> &dvals_x, const std::vector<double> &dvals_y){
            this->msh = msh;
            this->mu = mu;
            this->nv = msh.coords.nrows();
            this->nelems = msh.elems.nrows();
            this->ndofs = 3*nv+2*nelems;
            this->A = Eigen::SparseMatrix<double>(this->ndofs, this->ndofs);
            A.reserve(Eigen::VectorXi::Constant(this->ndofs,100));
            this->b = Eigen::VectorXd::Zero(this->ndofs);
            this->sol = Eigen::VectorXd::Zero(this->ndofs);
            this->u_x = std::vector<double>(this->nv,0);
            this->u_y = std::vector<double>(this->nv,0);
            this->p = std::vector<double>(this->nv,0);
            this->dvals_x = dvals_x;
            this->dvals_y = dvals_y;
            this->frhs_x = frhs_x;
            this->frhs_y = frhs_y;
            this->eps = 1e-3;
        }

        /**
         * @brief Assemble matrix using B1 space for velocity and P1 for pressure
         *
         *  | mu*grad(u1)*grad(v1)   0                      -p*div(v1) |
         *  | 0                      mu*grad(u2)*grad(v2)   -p*div(v2) |
         *  | -q*div(u1)             -q*div(u2)             -eps*p*q   |
         * 
         */
        void assemble_B1_P1();

        void assemble_3P1_P1();

        void apply_dbc();

        void solve();

        std::vector<double> u_x;
        std::vector<double> u_y;
        std::vector<double> p;
        std::vector<double> dvals_x;
        std::vector<double> dvals_y;
        std::vector<double> frhs_x;
        std::vector<double> frhs_y;
        Mesh2D msh;

    private:
        
        double mu;
        double eps;
        int nv;
        int nelems;
        int ndofs;
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        Eigen::VectorXd sol;
};


#endif // STOKES_HPP