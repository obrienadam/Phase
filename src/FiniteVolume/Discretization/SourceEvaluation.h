#ifndef SOURCE_EVALUATION_H
#define SOURCE_EVALUATION_H

#include "VectorFiniteVolumeField.h"

namespace fv
{

enum SourceEvaluationMethod{POTENTIAL};

VectorFiniteVolumeField source(VectorFiniteVolumeField field);
VectorFiniteVolumeField inverseWeightedSource(const ScalarFiniteVolumeField& w, const VectorFiniteVolumeField& field);

VectorFiniteVolumeField gravity(const ScalarFiniteVolumeField &rho, const Vector2D& g);

}

#endif
