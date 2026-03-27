# Math tests executable
add_executable(test_math MathTests.cpp)
target_link_libraries(test_math PRIVATE phase_math phase_system)
add_test(NAME MathTests COMMAND test_math)

# Geometry tests executable
add_executable(test_geometry GeometryTests.cpp)
target_link_libraries(test_geometry PRIVATE phase_2d_geometry phase_system)
add_test(NAME GeometryTests COMMAND test_geometry)
# Lid-driven cavity integration test
add_executable(test_lid_driven_cavity LidDrivenCavityTest.cpp)
add_test(NAME LidDrivenCavityTest COMMAND test_lid_driven_cavity)

# Backward-facing step integration test
add_executable(test_backward_step BackwardFacingStepTest.cpp)
add_test(NAME BackwardFacingStepTest COMMAND test_backward_step)