#include "TriMesh.hpp"

// constructor which triangulates set of coordinates
Mesh2D::Mesh2D(Matrix<double> &coords){

    int rows = coords.nrows();
    int cols = coords.ncols();
    DoubleMatrix xs = DoubleMatrix_create(rows, cols);

    for (int i = 0; i<rows; ++i){
        for (int j = 0; j<cols; ++j){
            xs.data[i*cols + j] = coords(i,j);
        }
    }

    Mesh msh = GeoMesh_Delaunay(&xs,1);

    this->coords = coords;

    this->elems = Matrix<int>(msh.nelems,3);
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<3; j++){
            this->elems(i,j) = msh.elems.data[3*i+j];
        }
    }

    this->sibhfs = Matrix<int>(msh.nelems,3);
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<3; j++){
            this->sibhfs(i,j) = msh.sibhfs.data[3*i+j];
        }
    }

    this->elems_on_boundary.resize(msh.nelems);
    for (int i = 0; i<msh.nelems; i++){
        this->elems_on_boundary[i] = msh.on_boundary[i];
    }

    Mesh_delete(&msh);
}

void Mesh2D::Triangulate(Matrix<double> &coords){
    int rows = coords.nrows();
    int cols = coords.ncols();
    DoubleMatrix xs = DoubleMatrix_create(rows, cols);

    for (int i = 0; i<rows; ++i){
        for (int j = 0; j<cols; ++j){
            xs.data[i*cols + j] = coords(i,j);
        }
    }

    Mesh msh = GeoMesh_Delaunay(&xs,1);

    this->coords = coords;

    this->elems = Matrix<int>(msh.nelems,3);
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<3; j++){
            this->elems(i,j) = msh.elems.data[3*i+j];
        }
    }

    this->sibhfs = Matrix<int>(msh.nelems,3);
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<3; j++){
            this->sibhfs(i,j) = msh.sibhfs.data[3*i+j];
        }
    }

    this->elems_on_boundary.resize(msh.nelems);
    for (int i = 0; i<msh.nelems; i++){
        this->elems_on_boundary[i] = msh.on_boundary[i];
    }

    Mesh_delete(&msh);
}


Mesh2D::Mesh2D(Matrix<double> &coords, Matrix<int> &segments){
    int rows = coords.nrows();
    int cols = coords.ncols();
    DoubleMatrix xs = DoubleMatrix_create(rows, cols);

    for (int i = 0; i<rows; ++i){ // copying coordinate array
        for (int j = 0; j<cols; ++j){
            xs.data[i*cols + j] = coords(i,j);
        }
    }

    rows = segments.nrows();
    cols = segments.ncols();
    IntMatrix segs = IntMatrix_create(rows, cols);

    for (int i = 0; i<rows; ++i){ // copying coordinate array
        for (int j = 0; j<cols; ++j){
            segs.data[i*cols + j] = segments(i,j);
        }
    }



    Mesh msh = GeoMesh_ConstrainedDelaunay(&segs, &xs);

    this->coords = coords;

    this->elems = Matrix<int>(msh.nelems,3);
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<3; j++){
            this->elems(i,j) = msh.elems.data[3*i+j];
        }
    }

    this->sibhfs = Matrix<int>(msh.nelems,3);
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<3; j++){
            this->sibhfs(i,j) = msh.sibhfs.data[3*i+j];
        }
    }

    this->elems_on_boundary.resize(msh.nelems);
    for (int i = 0; i<msh.nelems; i++){
        this->elems_on_boundary[i] = msh.on_boundary[i];
    }

    Mesh_delete(&msh);
}

void Mesh2D::Triangulate(Matrix<double> &coords, Matrix<int> &segments){
    int rows = coords.nrows();
    int cols = coords.ncols();
    DoubleMatrix xs = DoubleMatrix_create(rows, cols);

    for (int i = 0; i<rows; ++i){ // copying coordinate array
        for (int j = 0; j<cols; ++j){
            xs.data[i*cols + j] = coords(i,j);
        }
    }

    rows = segments.nrows();
    cols = segments.ncols();
    IntMatrix segs = IntMatrix_create(rows, cols);

    for (int i = 0; i<rows; ++i){ // copying coordinate array
        for (int j = 0; j<cols; ++j){
            segs.data[i*cols + j] = segments(i,j);
        }
    }



    Mesh msh = GeoMesh_ConstrainedDelaunay(&segs, &xs);

    this->coords = coords;

    this->elems = Matrix<int>(msh.nelems,3);
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<3; j++){
            this->elems(i,j) = msh.elems.data[3*i+j];
        }
    }

    this->sibhfs = Matrix<int>(msh.nelems,3);
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<3; j++){
            this->sibhfs(i,j) = msh.sibhfs.data[3*i+j];
        }
    }

    this->elems_on_boundary.resize(msh.nelems);
    for (int i = 0; i<msh.nelems; i++){
        this->elems_on_boundary[i] = msh.on_boundary[i];
    }

    Mesh_delete(&msh);
}


// printing function
void Mesh2D::print(){
    std::cout << "-------COORDINATES------" << std::endl;
    for (int i = 0; i<this->coords.nrows(); i++){
        std::cout << "(" << this->coords(i,0) << "," << this->coords(i,1) << ")" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "-------ELEMENTS------" << std::endl;
    for (int i = 0; i<this->coords.nrows(); i++){
        std::cout << this->elems(i,0) << "," << this->elems(i,1) << "," << this->elems(i,2) << " | " << this->elems_on_boundary[i] << std::endl;
    }
}