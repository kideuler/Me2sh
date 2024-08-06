cmake_minimum_required(VERSION 3.2)

project(
  Me2sh
  VERSION 0.0.1
  LANGUAGES CXX C)


# trimesh library
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Trimesh
    )

add_library(Trimesh 
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMesh_Linalg.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMesh2d.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMesh3d.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMeshPartition.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMeshSmooth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/kdTree.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/Matrix.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMesh.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/Trimesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/Trimesh.hpp
    )

 
add_executable(tm ${CMAKE_CURRENT_SOURCE_DIR}/Tests/TrimeshTest.cpp)

target_link_libraries(tm PRIVATE Trimesh)