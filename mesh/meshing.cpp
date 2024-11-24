#include <meshing.hpp>
#include <iostream>

void Me2sh_Mesh::generate(std::vector<int> PlaneTags, int MeshAlgo, int ElementType, int RecombinationAlgorithm, double h_target){
    int init = gmsh::isInitialized();
    if (init == 0){return;}

    gmsh::model::mesh::clear();

    // set the meshing algorithm
    gmsh::option::setNumber("Mesh.MeshSizeMin", h_target);
    gmsh::option::setNumber("Mesh.MeshSizeMax", 1.3*h_target);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.Algorithm", MeshAlgo);
    gmsh::option::setNumber("Mesh.RecombinationAlgorithm",RecombinationAlgorithm);

    // generate the mesh and recombine as needed
    if (ElementType == 2){
        gmsh::option::setNumber("Mesh.RecombineAll", 1);
    } else {
        gmsh::option::setNumber("Mesh.RecombineAll", 0);
    }
    gmsh::model::mesh::generate(2);

    // Get all entities
    gmsh::vectorpair entities;
    gmsh::model::getEntities(entities,2);
    int num_entities = entities.size();

    if (ElementType == 1){
        triMeshes.resize(num_entities);
    } else if (ElementType == 2){
        quadMeshes.resize(num_entities);
    }

    // Loop through all entities and get the 2D mesh for each entity
    std::vector<std::size_t> nodeTags;
    std::vector<double> Coords;
    std::vector<double> params;
    std::vector<int> map;

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t> > elementTags;
    std::vector<std::vector<std::size_t> > element_nodeTags;
    

    for (int i = 0; i < num_entities; i++){
        gmsh::model::mesh::rebuildNodeCache();
        gmsh::model::mesh::rebuildElementCache();
        int tag = entities[i].second;

        gmsh::model::mesh::getNodes(nodeTags, Coords, params, 2, tag, true, false);
        int num_nodes = nodeTags.size();

        gmsh::model::mesh::getElements(elementTypes, elementTags, element_nodeTags, 2, tag);

        std::size_t tritag = 0;
        std::size_t quadtag = 0;
        for (int j = 0; j < elementTypes.size(); j++){
            if (elementTypes[j] == 2){
                tritag = j;
            } else if (elementTypes[j] == 3){
                quadtag = j;
            }
        }

        
        if (ElementType == 1){
            triMeshes[i].second.resize(num_nodes);
            std::size_t nnodes;
            gmsh::model::mesh::getMaxNodeTag(nnodes);
            map.resize(nnodes+1);
            for (int j = 0; j < num_nodes; j++){
               triMeshes[i].second[j] = {Coords[3*j], Coords[3*j+1]};
               map[nodeTags[j]] = j;
            }

            int nelems = elementTags[tritag].size();
            triMeshes[i].first.resize(nelems);
            for (int j = 0; j < nelems; j++){
                triMeshes[i].first[j] = {map[element_nodeTags[tritag][3*j]],map[element_nodeTags[tritag][3*j+1]],map[element_nodeTags[tritag][3*j+2]]};
            }



        } else if (ElementType == 2){
            quadMeshes[i].second.resize(num_nodes);
            std::size_t nnodes;
            gmsh::model::mesh::getMaxNodeTag(nnodes);
            map.resize(nnodes+1);
            for (int j = 0; j < num_nodes; j++){
               quadMeshes[i].second[j] = {Coords[3*j], Coords[3*j+1]};
               map[nodeTags[j]] = j;
            }
            
            int nelems = elementTags[quadtag].size();
            quadMeshes[i].first.resize(nelems);
            for (int j = 0; j < nelems; j++){
                quadMeshes[i].first[j] = {map[element_nodeTags[quadtag][4*j]],map[element_nodeTags[quadtag][4*j+1]],map[element_nodeTags[quadtag][4*j+2]],map[element_nodeTags[quadtag][4*j+3]]};
            }

        }
    }

}