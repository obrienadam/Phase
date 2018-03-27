project(Phase)
cmake_minimum_required(VERSION 3.9.0)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -xHost")

    # Look for mkl by default
    set(BLA_VENDOR Intel10_64lp)
else ()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-reorder -Wno-sign-compare -Wno-switch -march=native")

    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "Requires at least gcc-4.9. You have gcc-${CMAKE_CXX_COMPILER_VERSION}.")
    endif ()

    # Look for OpenBLAS by default (includes CBLAS)
    set(BLA_VENDOR OpenBLAS)

    # Check for lapacke/cblas availability
    find_path(OPENBLAS_INCLUDE_DIRS NAMES "lapacke.h" "cblas.h" PATH_SUFFIXES "openblas")
    if (NOT OPENBLAS_INCLUDE_DIRS)
        message(FATAL_ERROR "lapacke.h and/or cblas.h headers not found.")
    endif ()

    include_directories(${OPENBLAS_INCLUDE_DIRS})
endif ()

# find required packages
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Boost 1.60 REQUIRED COMPONENTS filesystem system)
find_package(MPI REQUIRED)
find_package(Trilinos REQUIRED COMPONENTS Tpetra Belos MueLu Amesos2)
find_package(HDF5 REQUIRED)
find_package(OpenMP)

include_directories(${MPI_CXX_INCLUDE_PATH})

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Compiler configuration
set(CMAKE_CXX_STANDARD 11)

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

add_subdirectory(src)
add_subdirectory(utilities)


