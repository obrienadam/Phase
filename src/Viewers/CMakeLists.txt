set(HEADERS Viewer.h
            CgnsViewer.h
            IbViewer.h)

set(SOURCE  Viewer.cpp
            CgnsViewer.cpp
            IbViewer.cpp)

add_library(Viewers ${HEADERS} ${SOURCE})
target_link_libraries(Viewers boost_system boost_filesystem cgns)
