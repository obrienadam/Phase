#ifndef PHASE_FIXED_BOUNDARY_CONDITION_H
#define PHASE_FIXED_BOUNDARY_CONDITION_H

#include "BoundaryCondition.h"

template <class T> class FixedBoundaryCondition : public BoundaryCondition<T> {
public:
  Type type() const override { return FIXED; }

  const T &refVal() const { return _refVal; }

protected:
  T _refVal;
};

#endif
