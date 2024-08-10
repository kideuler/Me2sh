#ifndef FEM_HPP
#define FEM_HPP

#include <Eigen/Sparse>

// only dirichlet boundary conditions for now
class FemPoisson {

    public:
        FemPoisson(){};

        init(Mesh2D &msh, const std::vector<double> &frhs, double kappa, const std::vector<double> &dvals);
        init(Mesh2D &msh, const std::vector<double> &frhs, double kappa, double dval);
        init(Mesh2D &msh, double frhs, double kappa, const std::vector<double> &dvals);
        init(Mesh2D &msh, double frhs, double kappa, double dval);

        void assemble();

        void solve();
    
};


class FemEikonal {

    public:
        FemEikonal(){};

        init(Mesh2D &msh, double alpha);

        void assemble();

        void solve();
    
};

#endif