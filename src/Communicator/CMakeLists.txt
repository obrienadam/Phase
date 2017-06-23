set(HEADERS Communicator.h)
set(SOURCES Communicator.cpp)

add_library(Communicator ${HEADERS} ${SOURCES})
target_link_libraries(Communicator ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
