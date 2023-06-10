#ifndef PHASE_SLIP_BOUNDARY_CONDITION_H
#define PHASE_SLIP_BOUNDARY_CONDITION_H

#include "Geometry/Vector3D.h"

#include "BoundaryCondition.h"

class SlipBoundaryCondition : public BoundaryCondition<Vector3D> {
public:
  Type type() const override { return SLIP; }
};

#endif
