cmake_minimum_required(VERSION 3.20)
project(openFEM VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

find_package(Eigen3 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(fmt REQUIRED)
find_package(Threads REQUIRED)

# OpenCASCADE
find_library(OCC_KERNEL   TKernel    REQUIRED)
find_library(OCC_MATH     TKMath     REQUIRED)
find_library(OCC_BREP     TKBRep     REQUIRED)
find_library(OCC_STEP     TKSTEP     REQUIRED)
find_library(OCC_TOPOALGO TKTopAlgo  REQUIRED)
find_library(OCC_STEPBASE TKSTEPBase REQUIRED)
find_library(OCC_STEPATTR TKSTEPAttr REQUIRED)
find_library(OCC_XCAF     TKXCAF     REQUIRED)
find_library(OCC_XDESTEP  TKXDESTEP  REQUIRED)
find_library(OCC_GEOMBASE TKGeomBase REQUIRED)
find_library(OCC_G3D      TKG3d      REQUIRED)
find_library(OCC_XSBASE   TKXSBase   REQUIRED)

set(OCC_LIBS
    ${OCC_KERNEL} ${OCC_MATH}     ${OCC_BREP}
    ${OCC_STEP}   ${OCC_TOPOALGO} ${OCC_STEPBASE}
    ${OCC_STEPATTR} ${OCC_XCAF}   ${OCC_XDESTEP}
    ${OCC_GEOMBASE} ${OCC_G3D}    ${OCC_XSBASE}
)

find_path(OCC_INCLUDE "Standard.hxx"
          PATHS /usr/include/opencascade
                /usr/local/include/opencascade
          REQUIRED)

# gmsh
find_library(GMSH_LIB gmsh REQUIRED)
find_path(GMSH_INCLUDE gmsh.h
          PATHS /usr/include /usr/local/include
          REQUIRED)

file(GLOB SOURCES
    "src/cli/*.cpp"
    "src/elements/*.cpp"
    "src/fem/*.cpp"
    "src/geometry/*.cpp"
    "src/io/*.cpp"
    "src/mesh/*.cpp"
    "src/main.cpp"
)

add_executable(openFEM ${SOURCES})

target_include_directories(openFEM PRIVATE
    include
    ${OCC_INCLUDE}
    ${GMSH_INCLUDE}
)

target_link_libraries(openFEM PRIVATE
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    fmt::fmt
    Threads::Threads
    ${OCC_LIBS}
    ${GMSH_LIB}
)

add_executable(test_beam_frame
    src/tests/test_beam_frame.cpp
    src/elements/BeamElement.cpp
    src/elements/ShellElement.cpp
    src/elements/TetraElement.cpp
    src/fem/Assembler.cpp
    src/fem/Solver.cpp
)

target_include_directories(test_beam_frame PRIVATE
    include
    ${OCC_INCLUDE}
    ${GMSH_INCLUDE}
)

target_link_libraries(test_beam_frame PRIVATE
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    fmt::fmt
    Threads::Threads
    ${OCC_LIBS}
    ${GMSH_LIB}
)