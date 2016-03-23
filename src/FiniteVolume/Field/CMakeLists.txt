set(HEADERS Field.h
            IntegerField.h
            ScalarField.h
            VectorField.h
            PointField.h
            ScalarFiniteVolumeField.h
            VectorFiniteVolumeField.h)

set(SOURCES Field.tpp
            ScalarFiniteVolumeField.cpp
            VectorFiniteVolumeField.cpp)

add_library(Field ${HEADERS} ${SOURCES})
target_link_libraries(Field Geometry FiniteVolumeGrid2D)
