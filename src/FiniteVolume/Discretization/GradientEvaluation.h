#ifndef GRADIENT_EVALUATION_H
#define GRADIENT_EVALUATION_H

#include <functional>

#include "VectorFiniteVolumeField.h"

namespace fv
{

enum GradientEvaluationMethod{GREEN_GAUSS_CELL_CENTERED, GREEN_GAUSS_NODE_CENTERED, LEAST_SQUARES_CELL_CENTERED, FACE_TO_CELL};

//- Scalar gradients
VectorFiniteVolumeField computeGradient(GradientEvaluationMethod method, ScalarFiniteVolumeField& field, bool useCurrentFacesValues = false);
void computeGradient(GradientEvaluationMethod method, ScalarFiniteVolumeField& field, VectorFiniteVolumeField& gradField, bool useCurrentFaceValues = false);

//- Vector gradients (jacobians)

}

#endif
