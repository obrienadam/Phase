set(HEADERS Viewer.h)

set(SOURCE Viewer.cpp)

add_library(Viewers ${HEADERS} ${SOURCE})
target_link_libraries(Viewers cgns)
