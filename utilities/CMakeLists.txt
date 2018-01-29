add_executable(phaseReconstructSolution phaseReconstructSolution.cpp)
target_link_libraries(phaseReconstructSolution System ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} cgns hdf5)

add_executable(phasePartitionMesh phasePartitionMesh.cpp)
target_link_libraries(phasePartitionMesh FiniteVolumeGrid2D ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
