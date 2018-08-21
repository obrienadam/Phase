set(HEADERS Input.h
        CommandLine.h
        Exception.h
        StaticVector.h
        Communicator.h
        Timer.h
        RunControl.h
        NotImplementedException.h
        CgnsFile.h
        SolverInterface.h
        PostProcessingInterface.h)

set(SOURCES Input.cpp
        CommandLine.cpp
        Exception.cpp
        StaticVector.tpp
        Communicator.cpp
        Timer.cpp
        RunControl.cpp
        CgnsFile.cpp
        PostProcessingInterface.cpp)

add_library(phase_system ${HEADERS} ${SOURCES})

target_link_libraries(phase_system
        ${Boost_PROGRAM_OPTIONS_LIBRARY}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${MPI_C_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        cgns)

install(TARGETS
        phase_system
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
