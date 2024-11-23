#ifndef MESHING_H
#define MESHING_H

#include "gmsh.h"
#include <array>
#include <utility>

// Primitive types
typedef std::array<double, 2> Point;
typedef std::array<int, 2> Segment;
typedef std::array<int, 3> Triangle;
typedef std::array<int, 4> Quadrilateral;

// Mesh Types
typedef std::pair<std::vector<Triangle>, std::vector<Point>> TriMesh;
typedef std::pair<std::vector<Quadrilateral>, std::vector<Point>> QuadMesh;

// Partitioned Mesh types
typedef std::vector<TriMesh> PartitionedTriMesh;
typedef std::vector<QuadMesh> PartitionedQuadMesh;


class Me2sh_Mesh {
    public:
        Me2sh_Mesh(){ };
        ~Me2sh_Mesh(){};

        void clear(){
            triMeshes.clear();
            quadMeshes.clear();
        }

        void generate(std::vector<int> PlaneTags, int MeshAlgo, int ElementType, int RecombinationAlgorithm, double h_target);

        std::vector<TriMesh> triMeshes;
        std::vector<QuadMesh> quadMeshes;
    
};

#endif