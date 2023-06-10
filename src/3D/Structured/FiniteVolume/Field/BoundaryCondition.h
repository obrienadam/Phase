#ifndef PHASE_BOUNDARY_CONDITION_H
#define PHASE_BOUNDARY_CONDITION_H

#include "Structured/StructuredGrid3D/BoundaryPatch.h"

template <class T> class BoundaryCondition {
public:
  enum Type { FIXED, ZERO_GRADIENT, SLIP };

  virtual Type type() const = 0;
};

#endif
