set(HEADERS Viewer.h
        CgnsViewer.h
        PostProcessing.h
        PostProcessingObject.h
        ImmersedBoundaryObjectContactLineTracker.h
        IbTracker.h
        ImmersedBoundaryObjectProbe.h)

set(SOURCE Viewer.cpp
        CgnsViewer.cpp
        PostProcessing.cpp
        PostProcessingObject.cpp
        ImmersedBoundaryObjectContactLineTracker.cpp
        IbTracker.cpp
        ImmersedBoundaryObjectProbe.cpp)

add_library(PostProcessing ${HEADERS} ${SOURCE})
target_link_libraries(PostProcessing FiniteVolumeGrid2D ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

install(TARGETS
        PostProcessing
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)