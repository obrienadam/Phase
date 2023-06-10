#ifndef PHASE_BOUNDARY_CONDITION_FACTORY_H
#define PHASE_BOUNDARY_CONDITION_FACTORY_H

#include <memory>

#include "BoundaryCondition.h"

template <class T> class BoundaryConditionFactory {
public:
  std::unique_ptr<BoundaryCondition<T>> create(const std::string bcType);
};

#include "BoundaryConditionFactory.tpp"

#endif
