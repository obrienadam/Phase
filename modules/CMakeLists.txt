add_executable(phasePoisson phasePoisson.cpp)
target_link_libraries(phasePoisson ${PHASE_LIBRARIES})

add_executable(phasePiso phasePiso.cpp)
target_link_libraries(phasePiso ${PHASE_LIBRARIES})

add_executable(phasePisoMultiphase phasePisoMultiphase.cpp)
target_link_libraries(phasePisoMultiphase ${PHASE_LIBRARIES})

# This target has been deprecated
#add_executable(phaseCoupled phaseCoupled.cpp)
#target_link_libraries(phaseCoupled ${PHASE_LIBRARIES})

add_executable(phaseFractionalStep phaseFractionalStep.cpp)
target_link_libraries(phaseFractionalStep ${PHASE_LIBRARIES})

add_executable(phaseFractionalStepQuadraticIbm phaseFractionalStepQuadraticIbm.cpp)
target_link_libraries(phaseFractionalStepQuadraticIbm ${PHASE_LIBRARIES})

add_executable(phaseFractionalStepMultiphase phaseFractionalStepMultiphase.cpp)
target_link_libraries(phaseFractionalStepMultiphase ${PHASE_LIBRARIES})

add_executable(phaseFractionalStepIncremental phaseFractionalStepIncremental.cpp)
target_link_libraries(phaseFractionalStepIncremental ${PHASE_LIBRARIES})

add_executable(phaseFractionalStepIncrementalMultiphase phaseFractionalStepIncrementalMultiphase.cpp)
target_link_libraries(phaseFractionalStepIncrementalMultiphase ${PHASE_LIBRARIES})
