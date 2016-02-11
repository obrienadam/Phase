project(Phase)
cmake_minimum_required(VERSION 2.8)

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wno-reorder -Wno-sign-compare -Wno-switch -fopenmp -std=c++0x -Og")
set(CMAKE_CXX_FLAGS_RELEASE "-Wno-reorder -Wno-sign-compare -Wno-switch -fopenmp -std=c++0x -Ofast")

# Includes
set(INCLUDE_DIR src/Types
                src/System
                src/Geometry
                src/FiniteVolumeGrid2D
                src/FiniteVolume
                src/FiniteVolume/Field
                src/FiniteVolume/Equation
                src/Solvers
                src/Viewers)

include_directories(${INCLUDE_DIR})

add_subdirectory(src)
add_subdirectory(modules)
