add_executable(phaseReconstructSolution phaseReconstructSolution.cpp)
target_link_libraries(phaseReconstructSolution System ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} cgns ${HDF5_LIBRARIES})

add_executable(phasePartitionMesh phasePartitionMesh.cpp)
target_link_libraries(phasePartitionMesh FiniteVolumeGrid2D ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

install(TARGETS
        phaseReconstructSolution
        phasePartitionMesh
        RUNTIME DESTINATION bin)

