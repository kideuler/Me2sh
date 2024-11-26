cmake_minimum_required(VERSION 3.30.0)

set(CMAKE_OSX_ARCHITECTURES "arm64")
set(CMAKE_CXX_STANDARD 17)

project(
  Me2sh
  VERSION 0.0.5
  LANGUAGES CXX C)


# Define the version number as a preprocessor definition
add_definitions(-DPROJECT_VERSION="${PROJECT_VERSION}")

# add RPATH to installation
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # set build type
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)
endif()

# Enable sanitizers if build type is Debug
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
message(STATUS "Sanitizers enabled for Debug build")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address,undefined")
endif()


include(CTest)

# Find qt6
find_package(Qt6
    REQUIRED COMPONENTS Core Gui Widgets
)

# finding python
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

# Include Python headers and libraries
include_directories(${Python3_INCLUDE_DIRS})
link_directories(${Python3_LIBRARY_DIRS})

# finding gmsh
find_library(GMSH_LIB gmsh)
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()

find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()
message(STATUS "found gmsh at GMSH_INC: ${GMSH_INC}, GMSH_LIB: ${GMSH_LIB}")
include_directories(${GMSH_INC})


# find mfem
find_library(MFEM_LIBRARIES mfem)
if(NOT MFEM_LIBRARIES)
  message(FATAL_ERROR "Could not find libmfem")
endif()
find_path(MFEM_INCLUDE_DIRS mfem.hpp)
if(NOT MFEM_INCLUDE_DIRS)
  message(FATAL_ERROR "Could not find mfem.hpp")
endif()
message(STATUS "found mfem at MFEM_INCLUDE_DIRS: ${MFEM_INCLUDE_DIRS}, MFEM_LIBRARIES: ${MFEM_LIBRARIES}")
# find_library(mfem REQUIRED NAMES MFEM)
include_directories(${MFEM_INCLUDE_DIRS})

# find suitesparse
find_library(SUITESPARSE_LIBRARIES umfpack)
if(NOT SUITESPARSE_LIBRARIES)
  message(FATAL_ERROR "Could not find libsuitesparse")
endif()
find_path(SUITESPARSE_INCLUDE_DIRS suitesparse/umfpack.h)
if(NOT SUITESPARSE_INCLUDE_DIRS)
  message(FATAL_ERROR "Could not find umfpack.h")
endif()
message(STATUS "found suitesparse at SUITESPARSE_INCLUDE_DIRS: ${SUITESPARSE_INCLUDE_DIRS}, SUITESPARSE_LIBRARIES: ${SUITESPARSE_LIBRARIES}")
include_directories(${SUITESPARSE_INCLUDE_DIRS}/suitesparse)

# geometry library
file( GLOB GEOMETRY_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/geom/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/geom/*.hpp )
add_library(geometry SHARED ${GEOMETRY_SRCS})
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geom)
target_link_libraries(geometry PUBLIC ${GMSH_LIB})

# meshing library
file( GLOB MESHING_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/mesh/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/mesh/*.hpp )
add_library(meshing SHARED ${MESHING_SRCS})
target_include_directories(meshing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mesh)
target_link_libraries(meshing PUBLIC geometry)

# adding tests
add_subdirectory(test)

set(CMAKE_AUTOMOC ON) 
set(CMAKE_AUTORCC ON) 
set(CMAKE_AUTOUIC ON)

# Glob all gui files together
file( GLOB GUI_SRCS 
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/*.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/*.hpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/geom/*.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/geom/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/mesh/*.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/mesh/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/sim/*.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/sim/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/etc/*.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/etc/*.hpp )

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/resources/me2sh_hummingbird_1.icns PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

# Add the executable
qt_add_executable(me2sh
    ${GUI_SRCS}
    ${CMAKE_CURRENT_SOURCE_DIR}/resources/me2sh_hummingbird_1.icns
)

set_target_properties(me2sh PROPERTIES
    WIN32_EXECUTABLE TRUE
    MACOSX_BUNDLE TRUE
    MACOSX_BUNDLE_ICON_FILE me2sh_hummingbird_1.icns
)


target_include_directories(me2sh PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/geom
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh
    ${CMAKE_CURRENT_SOURCE_DIR}/gui
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/geom
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/mesh
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/sim
    ${CMAKE_CURRENT_SOURCE_DIR}/gui/etc
    ${Python3_INCLUDE_DIRS}
)

target_link_libraries(me2sh PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    geometry
    meshing
    ${Python3_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${MFEM_LIBRARIES}
)


include(GNUInstallDirs)


# Installation rules
install(TARGETS me2sh geometry meshing
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    BUNDLE DESTINATION bin
    PUBLIC_HEADER DESTINATION include
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/resources DESTINATION share/me2sh)