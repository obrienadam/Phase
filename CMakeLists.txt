project(Phase)
cmake_minimum_required(VERSION 2.8)

# Output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#Compile options
option(PHASE_USE_EIGEN "Use the Eigen3 linear solvers library" ON)
option(PHASE_USE_MPI "Compile with MPI support" ON)
option(PHASE_USE_HYPRE "Use the HYPRE linear solvers library" ON)

# Compiler flag setup
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wno-reorder -Wno-sign-compare -Wno-switch -fopenmp -std=c++11 -frounding-math -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-Wno-reorder -Wno-sign-compare -Wno-switch -fopenmp -std=c++11 -frounding-math -O3")

elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -g -qopenmp -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -qopenmp -O3")
endif()

# find required packages
find_package(Boost REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(PHASE_USE_MPI)
    find_package(MPI REQUIRED)
    set(CMAKE_CXX_COMPILER "mpic++")
    set(CMAKE_C_COMPILER "mpicc")
    add_definitions(-DPHASE_USE_MPI)
endif()

# Includes
set(INCLUDE_DIR src/Types
    src/System
    src/Math
    src/Geometry
    src/FiniteVolumeGrid2D
    src/FiniteVolumeGrid2D/Node
    src/FiniteVolumeGrid2D/Cell
    src/FiniteVolumeGrid2D/Face
    src/FiniteVolumeGrid2D/Link
    src/FiniteVolume
    src/FiniteVolume/Field
    src/FiniteVolume/Equation
    src/FiniteVolume/Discretization
    src/FiniteVolume/ImmersedBoundary
    src/Solvers
    src/Communicator
    src/Solvers/Multiphase
    src/Solvers/Coupled
    src/Viewers
    src/Utilities)

include_directories(${INCLUDE_DIR})

message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
message(STATUS "Languages: ${languages}")
message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CXX compiler command: ${CMAKE_CXX_COMPILER}")
message(STATUS "C compiler: ${CMAKE_C_COMPILER_ID}")
message(STATUS "C compiler command: ${CMAKE_C_COMPILER}")
message(STATUS "C compiler flags: ${CMAKE_C_FLAGS}")
message(STATUS "CXX compiler flags: ${CMAKE_CXX_FLAGS}")

add_subdirectory(src)
add_subdirectory(modules)
