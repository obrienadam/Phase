# Build Matrix test unit
add_executable(MatrixTest MatrixTest.cpp)
target_link_libraries(MatrixTest Equation ImmersedBoundary)

# Build Polygon test unit
add_executable(PolygonTest PolygonTest.cpp)
target_link_libraries(PolygonTest Geometry)


# Build range search test unit
add_executable(RangeSearchTest RangeSearchTest.cpp)
target_link_libraries(RangeSearchTest Geometry FiniteVolumeGrid2D)

# Build immersed boundary test unit
add_executable(ImmersedBoundaryTest ImmersedBoundaryTest.cpp)
target_link_libraries(ImmersedBoundaryTest ImmersedBoundary)
