#ifndef BEZIER_HPP
#define BEZIER_HPP

#include "TriMesh.hpp"

class Bezier2D {

    public:
        /**
         * @brief Default constructor for a new Bezier 2 D object
         * 
         */
        Bezier2D(){};
        
        /**
         * @brief Initialize the spline with control points and degree
         * 
         * @param coords coordinates of control points
         */
        void init(Matrix<double> coords);

        /**
         * @brief evaluate the spline coordinates at parameter t
         * 
         * @param t parameter value
         * @return std::array<double,2> 
         */
        std::array<double,2> eval(double t);

        /**
         * @brief evaluate the spline normal at parameter t
         * 
         * @param t parameter value
         * @return std::array<double,2> 
         */
        std::array<double,2> normal(double t);

        /**
         * @brief Create coordinates and segemnts from the spline
         * 
         * @param h_target 
         * @param coords 
         * @param segments 
         */
        void create_segments(double h_target, Matrix<double> &coords, Matrix<int> &segments);

        /**
         * @brief reset the spline
         * 
         */
        void reset(){
            control_points = Matrix<double>(0,0);
            nv = 0;
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
        int nv;
        double arclength;
};

#endif // BEZIER_HPP