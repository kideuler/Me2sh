#ifndef FEM_HPP
#define FEM_HPP

#include <Eigen/Sparse>
#include "TriMesh.hpp"
#include <math.h>

// only dirichlet boundary conditions for now
class FemPoisson {

    public:
        /**
         * @brief Construct a new Fem Poisson object
         * 
         */
        FemPoisson(){};

        // public solution
        std::vector<double> u;

        /**
         * @brief Initialize the Poisson problem
         * 
         * @param msh mesh
         * @param frhs right hand side
         * @param kappa kappa
         * @param dvals dircihlet values
         */
        void init(Mesh2D &msh, const std::vector<double> &frhs, double kappa, const std::vector<double> &dvals){
            this->msh = msh;
            this->kappa = kappa;
            int nv = msh.coords.nrows();
            this->A = Eigen::SparseMatrix<double>(nv,nv);
            A.reserve(Eigen::VectorXi::Constant(nv,10));
            this->b = Eigen::VectorXd::Zero(nv);
            this->sol = Eigen::VectorXd::Zero(nv-msh.bndnodes.size());
            this->u = std::vector<double>(nv,0);
            this->dvals = dvals;
            this->frhs = frhs;
        }
        void init(Mesh2D &msh, double frhs, double kappa, double dval){
            this->msh = msh;
            this->kappa = kappa;
            int nv = msh.coords.nrows();
            this->A = Eigen::SparseMatrix<double>(nv,nv);
            A.reserve(Eigen::VectorXi::Constant(nv,10));
            this->b = Eigen::VectorXd::Zero(nv);
            this->sol = Eigen::VectorXd::Zero(nv-msh.bndnodes.size());
            this->u = std::vector<double>(nv,0);
            this->dvals = std::vector<double>(nv,dval);
            this->frhs = std::vector<double>(nv,frhs);
        }

        /**
         * @brief Assemble the Poisson problem
         * 
         */
        void assemble();

        /**
         * @brief Apply Dirichlet boundary conditions
         * 
         */
        void apply_dbc();

        /**
         * @brief Solve the Poisson problem
         * 
         */
        void solve();
        
        /**
         * @brief Reset the Poisson problem
         * 
         */
        void reset();
    
    private:
        Mesh2D msh;
        double kappa;
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        Eigen::VectorXd sol;
        std::vector<double> dvals;
        std::vector<double> frhs;
};


class FemEikonal {

    public:
        FemEikonal(){};

        // public solution
        std::vector<double> u;
        Eigen::VectorXd sol;

        void init(Mesh2D &msh, double alpha){
            this->msh = msh;
            this->alpha = alpha;
            int nv = msh.coords.nrows();
            this->A = Eigen::SparseMatrix<double>(nv,nv);
            A.reserve(Eigen::VectorXi::Constant(nv,10));
            this->b = Eigen::VectorXd::Zero(nv);
            this->sol = Eigen::VectorXd::Zero(nv-msh.bndnodes.size());
            this->u = std::vector<double>(nv,0);
        }

        void assemble();

        void apply_dbc();

        void solve();

        void reset();

    private:
        Mesh2D msh;
        double alpha;
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        
    
};

#endif