add_executable(phase-reconstruct-solution PhaseReconstructSolution.cpp)
target_link_libraries(phase-reconstruct-solution
    phase_system
    ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} cgns ${HDF5_LIBRARIES})

add_executable(phase-reconstruct-compact-solution PhaseReconstructCompactSolution.cpp)
target_link_libraries(phase-reconstruct-compact-solution
    phase_system
    ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} cgns ${HDF5_LIBRARIES})

install(TARGETS
    phase-reconstruct-solution
    phase-reconstruct-compact-solution
    RUNTIME DESTINATION bin)
