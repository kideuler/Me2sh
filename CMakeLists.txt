cmake_minimum_required(VERSION 3.30.0)

set(CMAKE_OSX_ARCHITECTURES "arm64")
set(CMAKE_CXX_STANDARD 17)

project(
  Me2sh
  VERSION 0.0.3
  LANGUAGES CXX C)




include(CTest)

# set build type

# Find all necessary libraries
find_package(Qt6
    REQUIRED COMPONENTS Core Gui Widgets
)
find_package(Gmsh REQUIRED)


# geometry library
file( GLOB GEOMETRY_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/geom/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/geom/*.hpp )
add_library(geometry SHARED ${GEOMETRY_SRCS})
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/geom)
target_link_libraries(geometry PUBLIC Gmsh::Gmsh)

# adding tests
add_subdirectory(test)

qt_standard_project_setup()

# Glob all gui files together
file( GLOB GUI_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/gui/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/gui/*.hpp )

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/resources/appicon.icns PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

# Add the executable
qt_add_executable(me2sh
    ${GUI_SRCS}
    ${CMAKE_CURRENT_SOURCE_DIR}/resources/appicon.icns
)

set_target_properties(me2sh PROPERTIES
    WIN32_EXECUTABLE TRUE
    MACOSX_BUNDLE TRUE
    MACOSX_BUNDLE_ICON_FILE appicon.icns
)


target_include_directories(me2sh PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/geom
    ${CMAKE_CURRENT_SOURCE_DIR}/gui
)

target_link_libraries(me2sh PRIVATE
    Qt6::Core
    Qt6::Gui
    Qt6::Widgets
    geometry
)


