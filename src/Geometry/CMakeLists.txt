set(HEADERS Geometry.h
            Vector2D.h
            Point2D.h
            Shape2D.h
            Circle.h
            Polygon.h
            Triangle.h
            BoundingBox.h)

set(SOURCES Vector2D.cpp
            Shape2D.cpp
            Circle.cpp
            Polygon.cpp
            Triangle.cpp
            BoundingBox.cpp)

add_library(Geometry ${HEADERS} ${SOURCES})
target_link_libraries(Geometry System CGAL gmp)
