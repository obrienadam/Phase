# Build Matrix test unit
add_executable(MatrixTest MatrixTest.cpp)
target_link_libraries(MatrixTest Equation ImmersedBoundary)

# Build Polygon test unit
add_executable(PolygonTest PolygonTest.cpp)
target_link_libraries(PolygonTest Geometry)
