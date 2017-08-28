set(HEADERS ContactLineTracker.h
        ForceIntegrator.h
        IbTracker.h)

set(SOURCE ContactLineTracker.cpp
        ForceIntegrator.cpp
        IbTracker.cpp)

add_library(PostProcessing ${HEADERS} ${SOURCE})