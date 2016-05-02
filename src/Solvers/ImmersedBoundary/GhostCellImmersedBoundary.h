#ifndef GHOST_CELL_IMMERSED_BOUNDARY_H
#define GHOST_CELL_IMMERSED_BOUNDARY_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"

namespace gc {

Equation<ScalarFiniteVolumeField> ib(const ImmersedBoundaryObject& ibObj, ScalarFiniteVolumeField& field);
Equation<VectorFiniteVolumeField> ib(const ImmersedBoundaryObject &ibObj, VectorFiniteVolumeField &field);

}

#endif
