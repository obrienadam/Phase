add_executable(phasePoisson phasePoisson.cpp)
target_link_libraries(phasePoisson Solvers)

#add_executable(phasePiso phasePiso.cpp)
#target_link_libraries(phasePiso Solvers)

#add_executable(phasePisoMultiphase phasePisoMultiphase.cpp)
#target_link_libraries(phasePisoMultiphase Solvers)

# This target has been deprecated
#add_executable(phaseCoupled phaseCoupled.cpp)
#target_link_libraries(phaseCoupled Solvers)

add_executable(phase phase.cpp)
target_link_libraries(phase Solvers)

add_executable(phaseFractionalStep phaseFractionalStep.cpp)
target_link_libraries(phaseFractionalStep Solvers)

add_executable(phaseFractionalStepAxisymmetric phaseFractionalStepAxisymmetric.cpp)
target_link_libraries(phaseFractionalStepAxisymmetric Solvers)

add_executable(phaseFractionalStepQuadraticIbm phaseFractionalStepQuadraticIbm.cpp)
target_link_libraries(phaseFractionalStepQuadraticIbm Solvers)

add_executable(phaseFractionalStepEulerLagrange phaseFractionalStepEulerLagrange.cpp)
target_link_libraries(phaseFractionalStepEulerLagrange Solvers)

add_executable(phaseFractionalStepMultiphase phaseFractionalStepMultiphase.cpp)
target_link_libraries(phaseFractionalStepMultiphase Solvers)

add_executable(phaseFractionalStepMultiphaseQuadraticIbm phaseFractionalStepMultiphaseQuadraticIbm.cpp)
target_link_libraries(phaseFractionalStepMultiphaseQuadraticIbm Solvers)

install(TARGETS
        phaseFractionalStep
        phaseFractionalStepMultiphase
        phaseFractionalStepQuadraticIbm
        phaseFractionalStepMultiphaseQuadraticIbm
        RUNTIME DESTINATION bin)
