add_executable(phasePoisson phasePoisson.cpp)
target_link_libraries(phasePoisson System FiniteVolumeGrid2D Viewers Solvers)

add_executable(phaseSimple phaseSimple.cpp)
target_link_libraries(phaseSimple System FiniteVolumeGrid2D Viewers Solvers)


