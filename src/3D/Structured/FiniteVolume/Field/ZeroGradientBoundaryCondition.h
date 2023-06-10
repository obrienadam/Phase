#ifndef PHASE_ZERO_GRADIENT_BOUNDARY_CONDITION_H
#define PHASE_ZERO_GRADIENT_BOUNDARY_CONDITION_H

#include "BoundaryCondition.h"

template <class T>
class ZeroGradientBoundaryCondition : public BoundaryCondition<T> {
public:
  Type type() const override { return ZERO_GRADIENT; }
};

#endif
