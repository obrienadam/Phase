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

add_executable(phaseFractionalStepSimple phaseFractionalStepSimple.cpp)
target_link_libraries(phaseFractionalStepSimple ${PHASE_LIBRARIES})

add_executable(phaseFractionalStepSimpleMultiphase phaseFractionalStepSimpleMultiphase.cpp)
target_link_libraries(phaseFractionalStepSimpleMultiphase ${PHASE_LIBRARIES})

add_executable(phaseFractionalStepMultiphase phaseFractionalStepMultiphase.cpp)
target_link_libraries(phaseFractionalStepMultiphase ${PHASE_LIBRARIES})
