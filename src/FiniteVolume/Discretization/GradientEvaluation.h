#ifndef GRADIENT_EVALUATION_H
#define GRADIENT_EVALUATION_H

#include <functional>

#include "VectorFiniteVolumeField.h"
#include "TensorFiniteVolumeField.h"

namespace fv
{

enum GradientEvaluationMethod{GREEN_GAUSS_CELL_CENTERED, GREEN_GAUSS_NODE_CENTERED, LEAST_SQUARES_CELL_CENTERED, FACE_TO_CELL};

//- Scalar gradients
void computeGradient(GradientEvaluationMethod method, ScalarFiniteVolumeField& field, VectorFiniteVolumeField& gradField, bool useCurrentFaceValues = false);
void computeInverseWeightedGradient(const ScalarFiniteVolumeField& w, ScalarFiniteVolumeField& field, VectorFiniteVolumeField& gradField);

//- Vector gradients (jacobians)
void computeGradient(GradientEvaluationMethod method, VectorFiniteVolumeField &field, TensorFiniteVolumeField& gradField, bool useCurrentFacesValues = false);

}

#endif
