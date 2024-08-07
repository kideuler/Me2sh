cmake_minimum_required(VERSION 3.5)

project(
  Me2sh
  VERSION 0.0.1
  LANGUAGES CXX C)

set(CMAKE_OSX_ARCHITECTURES "arm64")

# trimesh library
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh #triangular mesh generation
    ${CMAKE_CURRENT_SOURCE_DIR}/gui #gui using qt
    )

add_library(Trimesh 
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMesh_Linalg.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMesh2d.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMeshPartition.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMeshSmooth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/kdTree.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/Matrix.c
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/GeoMesh.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/Trimesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh/Trimesh.hpp
    )



# adding test for trimesh library
add_executable(TrimeshTest ${CMAKE_CURRENT_SOURCE_DIR}/Tests/TrimeshTest.cpp)
target_link_libraries(TrimeshTest PRIVATE Trimesh)




find_package(Qt6
    REQUIRED COMPONENTS Core Gui Widgets
    OPTIONAL_COMPONENTS PrintSupport
)

qt_standard_project_setup()

qt_add_executable(me2sh
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/gui_main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/mainwindow.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/mainwindow.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/drawmesh.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/drawmesh.hpp
)

set_target_properties(me2sh PROPERTIES
    WIN32_EXECUTABLE TRUE
    MACOSX_BUNDLE TRUE
)

target_link_libraries(me2sh PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Trimesh
)