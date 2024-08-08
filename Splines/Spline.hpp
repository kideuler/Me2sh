#ifndef SPLINE_HPP
#define SPLINE_HPP

#include "TriMesh.hpp"

class Spline2D {

    public:
        /**
         * @brief Default constructor for a new Spline 2 D object
         * 
         */
        Spline2D(){};
        
        /**
         * @brief Initialize the spline with control points and degree
         * 
         * @param coords coordinates of control points
         * @param degree degree of spline used
         */
        void init(Matrix<double> coords, int degree = 3);

        std::array<double,2> eval(double t);

    private: 
        Matrix<double> control_points;
        Matrix<double> xweights;
        Matrix<double> yweights;
        std::vector<double> params;
        int nv;
        int degree;

        void Cubic_spline();
};

#endif