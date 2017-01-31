project(Phase)
cmake_minimum_required(VERSION 2.8)

# Output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# find required packages
find_package(Boost REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

option(PHASE_USE_MPI "Compile with MPI support" ON)
if(PHASE_USE_MPI)
    find_package(MPI REQUIRED)
    set(CMAKE_CXX_COMPILER "mpic++")
    set(CMAKE_C_COMPILER "mpicc")
    add_definitions(-DPHASE_USE_MPI)
endif()

option(PHASE_USE_EIGEN "Use the Eigen3 linear solvers library" ON)
if(PHASE_USE_EIGEN)
    add_definitions(-DPHASE_USE_EIGEN)
endif()

option(PHASE_USE_HYPRE "Use the HYPRE linear solvers library" ON)
if(PHASE_USE_HYPRE)
    add_definitions(-DPHASE_USE_HYPRE)
endif()

option(PHASE_USE_PETSC "Use the HYPRE linear solvers library" ON)
if(PHASE_USE_PETSC)
    add_definitions(-DPHASE_USE_PETSC)
endif()

# Compiler flag setup
set(CMAKE_CXX_FLAGS "-std=c++11")

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wno-reorder -Wno-sign-compare -Wno-switch -fopenmp -O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-Wno-reorder -Wno-sign-compare -Wno-switch -O3")

elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -mkl -qopenmp -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-mkl -qopenmp -O3 -march=native -mtune=native -march=native")
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
add_subdirectory(utilities)
