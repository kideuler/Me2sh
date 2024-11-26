#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

#include "mfem.hpp"
#include "geometry.hpp"
#include "meshing.hpp"
#include <vector>
#include <array>

enum SimulationType {
    Poisson,
    Stokes,
    Helmholz,
    Diffuision,
    Wave_Equation,
    Navier_Stokes
};

class Me2sh_Simulation {
    public:
        // member functions
        Me2sh_Simulation(){};
        ~Me2sh_Simulation(){};

        void init(SimulationType sim_type, std::shared_ptr<Me2sh_Mesh> mesh);


        // member variables
        std::vector<double> scalar_field;
        std::vector<std::array<double,2>> vector_field;
};


#endif // __SIMULATION_HPP__