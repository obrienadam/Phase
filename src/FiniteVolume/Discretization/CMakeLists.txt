set(HEADERS CrankNicolson.h
            FaceInterpolation.h
            GradientEvaluation.h)

set(SOURCES CrankNicolson.cpp
            FaceInterpolation.tpp
            FaceInterpolation.cpp
            GradientEvaluation.cpp)

add_library(Discretization ${HEADERS} ${SOURCES})
target_link_libraries(Discretization Field Math)
