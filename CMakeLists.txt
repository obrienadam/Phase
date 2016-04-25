project(Phase)
cmake_minimum_required(VERSION 2.8)

# Output paths

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Compiler flags

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	
        set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wno-reorder -Wno-sign-compare -Wno-switch -fopenmp -std=c++11 -O0")
        set(CMAKE_CXX_FLAGS_RELEASE "-Wno-reorder -Wno-sign-compare -Wno-switch -fopenmp -std=c++11 -O3")
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})

elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")

	set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -g -qopenmp -O0")
        set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -qopenmp -O3")
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})

endif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

# Includes
set(INCLUDE_DIR src/Types
                src/System
                src/Geometry
                src/FiniteVolumeGrid2D
                src/FiniteVolume
                src/FiniteVolume/Field
                src/FiniteVolume/BoundaryCondition
                src/FiniteVolume/Equation
                src/Solvers
                src/Solvers/Multiphase
                src/Solvers/ImmersedBoundary
                src/Viewers)

include_directories(${INCLUDE_DIR})

message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
message(STATUS "Languages: ${languages}")
message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C compiler: ${CMAKE_C_COMPILER_ID}")
message(STATUS "CXX release flags: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "C release flags: ${CMAKE_C_FLAGS_RELEASE}")
message(STATUS "CXX debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "C debug flags: ${CMAKE_C_FLAGS_DEBUG}")

find_package(CGAL REQUIRED)

add_subdirectory(src)
add_subdirectory(modules)
add_subdirectory(test)
