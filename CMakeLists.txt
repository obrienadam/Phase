project(Phase)
cmake_minimum_required(VERSION 3.6.0)

# Output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Look for OpenBLAS by default (includes CBLAS)
set(BLA_VENDOR OpenBLAS)

# find required packages
find_package(Boost REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)
find_package(Trilinos REQUIRED)
find_package(HDF5 REQUIRED)

include_directories(${MPI_C_INCLUDE_PATH})
include_directories(${Trilinos_INCLUDE_DIRS})

# Check for lapacke availability
find_path(LAPACKE_INCLUDE_PATH "lapacke.h" PATHS "/usr/include/lapacke" "/usr/include/openblas")
if(NOT LAPACKE_INCLUDE_PATH)
    message(FATAL_ERROR "lapacke.h header not found.")
endif()

include_directories(${LAPACKE_INCLUDE_PATH})

# Compiler configuration
set(CMAKE_CXX_STANDARD 11)

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wno-reorder -Wno-sign-compare -Wno-switch -fopenmp -O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-Wno-reorder -Wno-sign-compare -Wno-switch -O3 -march=native -DNDEBUG")

    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "Requires at least gcc-4.9. You have gcc-${CMAKE_CXX_COMPILER_VERSION}.")
    endif()

elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -mkl -qopenmp -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-mkl -qopenmp -O3 -xHost -DNDEBUG")
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
        src/Viewers
        src/Utilities)

include_directories(${INCLUDE_DIR})

message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})
message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CXX compiler command: ${CMAKE_CXX_COMPILER}")
message(STATUS "C compiler: ${CMAKE_C_COMPILER_ID}")
message(STATUS "C compiler command: ${CMAKE_C_COMPILER}")
message(STATUS "C compiler flags: ${CMAKE_C_FLAGS}")
message(STATUS "CXX compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Boost include directory: " ${Boost_INCLUDE_DIRS})
message(STATUS "Boost library directory: " ${Boost_LIBRARY_DIRS})
message(STATUS "BLAS library: " ${BLAS_LIBRARIES})
message(STATUS "LAPACK library: " ${LAPACK_LIBRARIES})
message(STATUS "MPI include directory: " ${MPI_CXX_INCLUDE_PATH})
message(STATUS "MPI libraries: " ${MPI_C_LIBRARIES})
message(STATUS "Trilinos directory: " ${Trilinos_DIR})

set(PHASE_LIBRARIES Solvers FiniteVolumeGrid2D FiniteVolume Viewers Utilities Geometry Math System Communicator)

add_subdirectory(src)
add_subdirectory(modules)
add_subdirectory(utilities)
