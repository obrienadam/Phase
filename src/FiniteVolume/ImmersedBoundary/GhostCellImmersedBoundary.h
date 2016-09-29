#ifndef GHOST_CELL_IMMERSED_BOUNDARY_H
#define GHOST_CELL_IMMERSED_BOUNDARY_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"

namespace gc {

Equation<ScalarFiniteVolumeField> ib(const std::vector<ImmersedBoundaryObject>& ibObjs, ScalarFiniteVolumeField& field);

Equation<VectorFiniteVolumeField> ib(const std::vector<ImmersedBoundaryObject>& ibObjs, VectorFiniteVolumeField& field);
Equation<VectorFiniteVolumeField> ib(const std::vector<ImmersedBoundaryObject>& ibObjs, const ScalarFiniteVolumeField& rho, VectorFiniteVolumeField& field);

}

#endif
