cmake_minimum_required(VERSION 3.5)

project(
  Me2sh
  VERSION 0.0.3
  LANGUAGES CXX C)

include(CTest)

set(CMAKE_OSX_ARCHITECTURES "arm64")

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/Trimesh #triangular mesh generation
    ${CMAKE_CURRENT_SOURCE_DIR}/gui #gui using qt
    ${CMAKE_CURRENT_SOURCE_DIR}/Splines #splines for geometry
    )

# trimesh library
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


# finding eigen and adding spline library
find_package(Eigen3 REQUIRED)
add_library(Splines
    ${CMAKE_CURRENT_SOURCE_DIR}/Splines/Spline.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Splines/Spline.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Splines/Bezier.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Splines/Bezier.hpp
)
target_link_libraries(Splines Eigen3::Eigen)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/ConsoleOutput.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/ConsoleOutput.cpp
)

set_target_properties(me2sh PROPERTIES
    WIN32_EXECUTABLE TRUE
    MACOSX_BUNDLE TRUE
)

target_link_libraries(me2sh PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    Splines
    Trimesh
)


# adding test for trimesh library
add_executable(TrimeshTest ${CMAKE_CURRENT_SOURCE_DIR}/Tests/TrimeshTest.cpp)
target_link_libraries(TrimeshTest PRIVATE Splines Trimesh)
add_test(NAME TrimeshTest
         COMMAND TrimeshTest)

# adding test for spline library
add_executable(SplinemeshTest ${CMAKE_CURRENT_SOURCE_DIR}/Tests/SplinemeshTest.cpp)
target_link_libraries(SplinemeshTest PRIVATE Splines Trimesh)
add_test(NAME SplinemeshTest
         COMMAND SplinemeshTest)

# adding test for Bezier library
add_executable(BeziermeshTest ${CMAKE_CURRENT_SOURCE_DIR}/Tests/BeziermeshTest.cpp)
target_link_libraries(BeziermeshTest PRIVATE Splines Trimesh)
add_test(NAME BeziermeshTest
         COMMAND BeziermeshTest)