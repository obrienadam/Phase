#ifndef FACE_INTERPOLATION_H
#define FACE_INTERPOLATION_H

#include "FiniteVolumeField.h"

namespace fv
{

enum InterpolationMethod{UNWEIGHTED, INVERSE_VOLUME, INVERSE_DISTANCE, INVERSE_SQR_DISTANCE};

template<class T>
void interpolateFaces(InterpolationMethod method, FiniteVolumeField<T>& field);

template<class T>
void harmonicInterpolateFaces(InterpolationMethod method, FiniteVolumeField<T>& field);

}

#include "FaceInterpolation.tpp"

#endif
