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

void Mesh2D::Refine(double h_target){
    Mesh msh;
    msh.coords = DoubleMatrix_create(this->coords.nrows(), this->coords.ncols());
    msh.elems = IntMatrix_create(this->elems.nrows(), this->elems.ncols());
    msh.sibhfs = IntMatrix_create(this->sibhfs.nrows(), this->sibhfs.ncols());
    msh.nelems = this->elems.nrows();
    msh.on_boundary = (bool*) malloc(msh.nelems*sizeof(bool));
    msh.delete_elem = (bool*) malloc(msh.nelems*sizeof(bool));
    msh.bwork = (bool*) malloc(msh.nelems*sizeof(bool));
    msh.stack = (int*) malloc(msh.nelems*sizeof(int));
    msh.haskdTree = false;
    msh.haskdTree = false;
    msh.hasGraph = false;
    msh.hasPartition = false;

    for (int i = 0; i<this->coords.nrows(); i++){
        for (int j = 0; j<this->coords.ncols(); j++){
            msh.coords.data[i*this->coords.ncols() + j] = this->coords(i,j);
        }
    }

    for (int i = 0; i<this->elems.nrows(); i++){
        for (int j = 0; j<this->elems.ncols(); j++){
            msh.elems.data[i*this->elems.ncols() + j] = this->elems(i,j);
        }
    }

    for (int i = 0; i<this->sibhfs.nrows(); i++){
        for (int j = 0; j<this->sibhfs.ncols(); j++){
            msh.sibhfs.data[i*this->sibhfs.ncols() + j] = this->sibhfs(i,j);
        }
    }

    for (int i = 0; i<this->elems_on_boundary.size(); i++){
        msh.on_boundary[i] = this->elems_on_boundary[i];
    }

    if (h_target == 0){
        GeoMesh_DelaunayRefine(&msh, false, h_target, 2);
    } else {
        GeoMesh_DelaunayRefine(&msh, true, h_target, 2);
    }
    

    this->coords = Matrix<double>(msh.coords.nrows,2);
    for (int i = 0; i<msh.coords.nrows; i++){
        for (int j = 0; j<2; j++){
            this->coords(i,j) = msh.coords.data[2*i+j];
        }
    }

    this->elems = Matrix<int>(msh.elems.nrows,3);
    for (int i = 0; i<msh.elems.nrows; i++){
        for (int j = 0; j<3; j++){
            this->elems(i,j) = msh.elems.data[3*i+j];
        }
    }

    this->sibhfs = Matrix<int>(msh.sibhfs.nrows,3);
    for (int i = 0; i<msh.sibhfs.nrows; i++){
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

void Mesh2D::Smooth(int niters){
    Mesh msh;
    msh.coords = DoubleMatrix_create(this->coords.nrows(), this->coords.ncols());
    msh.elems = IntMatrix_create(this->elems.nrows(), this->elems.ncols());
    msh.sibhfs = IntMatrix_create(this->sibhfs.nrows(), this->sibhfs.ncols());
    msh.nelems = this->elems.nrows();
    msh.on_boundary = (bool*) malloc(msh.nelems*sizeof(bool));
    msh.delete_elem = (bool*) malloc(msh.nelems*sizeof(bool));
    msh.bwork = (bool*) malloc(msh.nelems*sizeof(bool));
    msh.stack = (int*) malloc(msh.nelems*sizeof(int));
    msh.haskdTree = false;
    msh.haskdTree = false;
    msh.hasGraph = false;
    msh.hasPartition = false;

    for (int i = 0; i<this->coords.nrows(); i++){
        for (int j = 0; j<this->coords.ncols(); j++){
            msh.coords.data[i*this->coords.ncols() + j] = this->coords(i,j);
        }
    }

    for (int i = 0; i<this->elems.nrows(); i++){
        for (int j = 0; j<this->elems.ncols(); j++){
            msh.elems.data[i*this->elems.ncols() + j] = this->elems(i,j);
        }
    }

    for (int i = 0; i<this->sibhfs.nrows(); i++){
        for (int j = 0; j<this->sibhfs.ncols(); j++){
            msh.sibhfs.data[i*this->sibhfs.ncols() + j] = this->sibhfs(i,j);
        }
    }

    for (int i = 0; i<this->elems_on_boundary.size(); i++){
        msh.on_boundary[i] = this->elems_on_boundary[i];
    }

    bool* bdy = Mesh_find_bdy_nodes(&msh);
    Mesh_smooth2d(&msh, bdy, niters);

    this->coords = Matrix<double>(msh.coords.nrows,2);
    for (int i = 0; i<msh.coords.nrows; i++){
        for (int j = 0; j<2; j++){
            this->coords(i,j) = msh.coords.data[2*i+j];
        }
    }

    this->elems = Matrix<int>(msh.elems.nrows,3);
    for (int i = 0; i<msh.elems.nrows; i++){
        for (int j = 0; j<3; j++){
            this->elems(i,j) = msh.elems.data[3*i+j];
        }
    }

    this->sibhfs = Matrix<int>(msh.sibhfs.nrows,3);
    for (int i = 0; i<msh.sibhfs.nrows; i++){
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

void Mesh2D::Compute_volume_length_metric(){
    quality.resize(this->elems.nrows());
    std::fill(quality.begin(), quality.end(), 0.0);
    double Max = -1e34;
    double Min = 1e34;
    for (int i = 0; i<elems.nrows(); i++){
        int i0 = elems(i,0);
        int i1 = elems(i,1);
        int i2 = elems(i,2);

        double x0 = coords(i0,0);
        double y0 = coords(i0,1);
        double x1 = coords(i1,0);
        double y1 = coords(i1,1);
        double x2 = coords(i2,0);
        double y2 = coords(i2,1);

        double area = 0.5*((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
        double l0 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
        double l1 = ((x0-x2)*(x0-x2) + (y0-y2)*(y0-y2));
        double l2 = ((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));

        quality[i] = (4*sqrt(3))*(area/(l0+l1+l2));
        if (quality[i] > Max){
            Max = quality[i];
        }
        if (quality[i] < Min){
            Min = quality[i];
        }
    }
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