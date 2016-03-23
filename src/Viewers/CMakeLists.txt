set(HEADERS Viewer.h
            TecplotViewer.h
            MatplotlibViewer.h)

set(SOURCE Viewer.cpp
           TecplotViewer.cpp
           MatplotlibViewer.cpp)

add_library(Viewers ${HEADERS} ${SOURCE})
target_link_libraries(Viewers python2.7)
