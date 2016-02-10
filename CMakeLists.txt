project(Phase)
cmake_minimum_required(VERSION 2.8)

# Compiler flags
# Compiler flags
set(CMAKE_CXX_FLAGS "-g -Wall -Wno-reorder -Wno-sign-compare -Wno-switch -fopenmp -std=c++0x -O0")
set(CMAKE_C_FLAGS "-Wall")

# Includes
set(INCLUDE_DIR src/Types
                src/InputOutput
                src/Geometry
                src/FiniteVolumeGrid2D
                src/FiniteVolume
                src/FiniteVolume/Field
                src/FiniteVolume/Equation
                src/Viewers)

include_directories(${INCLUDE_DIR})

add_subdirectory(src)
add_subdirectory(modules)
