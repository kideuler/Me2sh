#ifndef TRIMESH
#define TRIMESH

extern "C" {
#include "GeoMesh.h"
}

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>
#include <ctime>

/**
 * @brief Matrix class for handling 2D matrices
 * 
 * @tparam T data type of matrix
 */
template <typename T> class Matrix { 

    public:
        std::vector<T> arr;
        size_t rows;
        size_t cols;

        /**
         * @brief Construct a new Matrix object with default size 0,0
         * 
         */
        Matrix() : rows(0), cols(0) {}

        /**
         * @brief Construct a new Matrix object with specified size and set entries to 0
         * 
         * @param nrows number of rows
         * @param ncols number of columns
         */
        Matrix(size_t nrows, size_t ncols){
            rows = nrows;
            cols = ncols;
            arr.reserve(rows*cols);
            arr.resize(rows*cols);
            std::fill(arr.begin(), arr.end(), 0);
        }

        /**
         * @brief Fill matrix with random numbers between 0 and 1
         * 
         */
        void fill_rand(){
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0, 1);
            for (int i = 0; i<rows; ++i){
                for (int j = 0; j<cols; ++j){
                    arr[i*cols + j] = dis(gen);
                }
            }
        }

        /**
         * @brief Get number of rows
         * 
         * @return int number of rows
         */
        int nrows(){
            return rows;
        }

        /**
         * @brief Get number of columns
         * 
         * @return int number of columns
         */
        int ncols(){
            return cols;
        }

        /**
         * @brief Access matrix element at (i,j)
         * 
         * @param i row index
         * @param j column index
         * @return T& reference to element
         */
        T& operator()(size_t i, size_t j){
            return arr[i*cols + j];
        }

        /**
         * @brief Access matrix element at (i,j)
         * 
         * @param i row index
         * @param j column index
         * @return const T& reference to element
         */
        const T& operator()(size_t i, size_t j) const {
            return arr[i*cols + j];
        }
};

/**
 * @brief 2D mesh class
 * 
 */
class Mesh2D {

    public:
        Matrix<double> coords; // coordinates of mesh
        Matrix<int> elems; // elements of mesh
        Matrix<int> sibhfs; // sibling half-faces
        std::vector<bool> elems_on_boundary; // boolean array of elements on boundary
        std::vector<bool> nodes_on_boundary; // boolean array of nodes on boundary
        std::vector<double> quality; // quality of elements

        /**
         * @brief Construct a new Mesh 2D object with default size 0,0
         * 
         */
        Mesh2D(){
            coords = Matrix<double>(0,0);
            elems = Matrix<int>(0,0);
            sibhfs = Matrix<int>(0,0);
        }

        void reset(){
            coords = Matrix<double>(0,0);
            elems = Matrix<int>(0,0);
            sibhfs = Matrix<int>(0,0);
        }

        /**
         * @brief Construct a new Mesh 2D object from coordinates and triangulate them
         * 
         * @param coords coordinate array
         */
        Mesh2D(Matrix<double> &coords);

        /**
         * @brief triangulate a set of coordinates
         * 
         * @param coords coordinate array
         */
        void Triangulate(Matrix<double> &coords);

        /**
         * @brief Construct a new Mesh 2D object from coordinates and segments
         * 
         * @param coords  coordinate array
         * @param segments segment array
         */
        Mesh2D(Matrix<double> &coords, Matrix<int> &segments);

        /**
         * @brief constrained triangulation of a set of coordinates and segments
         * 
         * @param coords  coordinate array
         * @param segments segment array
         */
        void Triangulate(Matrix<double> &coords, Matrix<int> &segments);

        /**
         * @brief Refine mesh based on edge length
         * 
         * @param h_target target edge length
         */
        void Refine(double h_target);
        
        /**
         * @brief Smooth mesh with max number of iterations
         * 
         * @param niters 
         */
        void Smooth(int niters);

        /**
         * @brief compute volume length metric of the mesh
         * 
         */
        void Compute_volume_length_metric();

        /**
         * @brief print contents of the mesh
         * 
         */
        void print();
};


#endif