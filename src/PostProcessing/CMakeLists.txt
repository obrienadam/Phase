set(HEADERS PostProcessing.h
        ReynoldsNumber.h
        GhostCellImmersedBoundaryObjectContactLineTracker.h
        ContactLineTracker.h
        ForceIntegrator.h
        IbTracker.h)

set(SOURCE PostProcessing.cpp
        ReynoldsNumber.cpp
        GhostCellImmersedBoundaryObjectContactLineTracker.cpp
        ContactLineTracker.cpp
        ForceIntegrator.cpp
        IbTracker.cpp)

add_library(PostProcessing ${HEADERS} ${SOURCE})