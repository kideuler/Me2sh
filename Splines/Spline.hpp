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
        void init(const Matrix<double> &coords, int degree = 3);

        std::array<double,2> eval(double t);

        std::array<double,2> normal(double t);

        void create_segments(double h_target, Matrix<double> &coords, Matrix<int> &segments);

        void reset(){
            control_points = Matrix<double>(0,0);
            xweights = Matrix<double>(0,0);
            yweights = Matrix<double>(0,0);
            params = std::vector<double>(0);
            nv = 0;
            degree = 0;
            arclength = 0.0;
        }

        double get_arclength(){
            return arclength;
        }

        int npoints(){
            return nv;
        }

    private: 
        Matrix<double> control_points;
        Matrix<double> xweights;
        Matrix<double> yweights;
        std::vector<double> params;
        int nv;
        int degree;
        double arclength;

        void Cubic_spline();
};

#endif