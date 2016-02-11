set(HEADERS TecplotViewer.h
            MatplotlibViewer.h)

set(SOURCE TecplotViewer.cpp
           MatplotlibViewer.cpp)

add_library(Viewers ${HEADERS} ${SOURCE})
target_link_libraries(Viewers python2.7)
