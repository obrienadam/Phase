set(HEADERS Integrator.h
            ForceIntegrator.h
            VolumeIntegrator.h)

set(SOURCE Integrator.cpp
           ForceIntegrator.cpp
           VolumeIntegrator.cpp)

add_library(Utilities ${HEADERS} ${SOURCE})
target_link_libraries(Utilities Math)
