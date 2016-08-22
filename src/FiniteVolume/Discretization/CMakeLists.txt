set(HEADERS CrankNicolson.h
            AdamsBashforth.h
            FaceInterpolation.h
            GradientEvaluation.h)

set(SOURCES CrankNicolson.cpp
            AdamsBashforth.cpp
            FaceInterpolation.tpp
            GradientEvaluation.cpp)

add_library(Discretization ${HEADERS} ${SOURCES})
target_link_libraries(Discretization Field Math)
