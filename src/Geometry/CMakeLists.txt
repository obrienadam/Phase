set(HEADERS Vector2D.h
            Point2D.h
            Tensor2D.h
            Line2D.h
            LineSegment2D.h
            Ray2D.h
            Shape2D.h
            Circle.h
            Polygon.h
            BoundingBox.h
            Search.h
            Box.h)

set(SOURCES Vector2D.cpp
            Tensor2D.cpp
            Line2D.cpp
            LineSegment2D.cpp
            Ray2D.cpp
            Shape2D.cpp
            Circle.cpp
            Polygon.cpp
            BoundingBox.cpp
            Search.cpp
            Box.cpp)

add_library(Geometry ${HEADERS} ${SOURCES})
target_link_libraries(Geometry System Equation)
