set(HEADERS Vector2D.h
            Point2D.h
            Line2D.h
            Shape2D.h
            Circle.h
            Quadrilateral.h
            Polygon.h
            BoundingBox.h)

set(SOURCES Vector2D.cpp
            Line2D.cpp
            Shape2D.cpp
            Circle.cpp
            Quadrilateral.cpp
            Polygon.cpp
            BoundingBox.cpp)

add_library(Geometry ${HEADERS} ${SOURCES})
