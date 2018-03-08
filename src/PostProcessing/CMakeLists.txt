set(HEADERS Viewer.h
        CgnsViewer.h
        PostProcessing.h
        PostProcessingObject.h
        VolumeIntegrator.h
        ReynoldsNumber.h
        ImmersedBoundaryObjectContactLineTracker.h
        ForceIntegrator.h
        IbTracker.h)

set(SOURCE Viewer.cpp
        CgnsViewer.cpp
        PostProcessing.cpp
        PostProcessingObject.cpp
        VolumeIntegrator.cpp
        ReynoldsNumber.cpp
        ImmersedBoundaryObjectContactLineTracker.cpp
        ForceIntegrator.cpp
        IbTracker.cpp)

add_library(PostProcessing ${HEADERS} ${SOURCE})
target_link_libraries(PostProcessing FiniteVolumeGrid2D ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

install(TARGETS
        PostProcessing
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)