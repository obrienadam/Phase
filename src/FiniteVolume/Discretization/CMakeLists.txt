set(HEADERS CrankNicolson.h
            FaceInterpolation.h
            GradientEvaluation.h
            SourceEvaluation.h
            Cicsam.h
            Plic.h)

set(SOURCES CrankNicolson.cpp
            FaceInterpolation.tpp
            FaceInterpolation.cpp
            GradientEvaluation.cpp
            SourceEvaluation.cpp
            Cicsam.cpp
            Plic.cpp)

add_library(Discretization ${HEADERS} ${SOURCES})
target_link_libraries(Discretization Field Math)
