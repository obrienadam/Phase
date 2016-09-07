set(HEADERS Vector2D.h
            Point2D.h
            Tensor2D.h
            Line2D.h
            Shape2D.h
            Circle.h
            Polygon.h
            BoundingBox.h
            Search.h)

set(SOURCES Vector2D.cpp
            Tensor2D.cpp
            Line2D.cpp
            Shape2D.cpp
            Circle.cpp
            Polygon.cpp
            BoundingBox.cpp
            Search.tpp)

add_library(Geometry ${HEADERS} ${SOURCES})
target_link_libraries(Geometry System Equation)
