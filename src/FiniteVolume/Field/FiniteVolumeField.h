#ifndef FINITE_VOLUME_FIELD_H
#define FINITE_VOLUME_FIELD_H

#include "Field.h"
#include "FiniteVolumeGrid2D.h"

template <class T>
class FiniteVolumeField : public Field<T>
{
public:

    enum FaceInterpolationMethod{VOLUME_WEIGHTED, DISTANCE_WEIGHTED};

    FiniteVolumeField(const FiniteVolumeGrid2D& grid, const T& initialValue = T(), const std::string& name = "N/A");
    FiniteVolumeField(const FiniteVolumeGrid2D& grid, const std::string &name);

    const FiniteVolumeGrid2D& grid() const { return grid_; }

protected:

    const FiniteVolumeGrid2D& grid_;
    Field<T> facesI_, facesJ_;
    FaceInterpolationMethod faceInterpolationMethod_;

};

#include "FiniteVolumeField.tpp"

#endif
