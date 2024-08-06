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

template <typename T> class Matrix {
    private:
        std::vector<T> arr;
        size_t rows;
        size_t cols;

    public:

        // default constructor
        Matrix() : rows(0), cols(0) {}

        // initializing matrix
        Matrix(size_t nrows, size_t ncols){
            rows = nrows;
            cols = ncols;
            arr.reserve(rows*cols);
            arr.resize(rows*cols);
            std::fill(arr.begin(), arr.end(), 0);
        }

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

        int nrows(){
            return rows;
        }
        int ncols(){
            return cols;
        }

        // overloading access operator
        T& operator()(size_t i, size_t j){
            return arr[i*cols + j];
        }

        const T& operator()(size_t i, size_t j) const {
            return arr[i*cols + j];
        }
};

class Mesh2D {

    public:
        Matrix<double> coords;
        Matrix<int> elems;
        Matrix<int> sibhfs;

        Mesh2D(Matrix<double> &coords);
};


#endif