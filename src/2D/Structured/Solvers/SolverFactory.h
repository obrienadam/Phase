#ifndef PHASE_SOLVER_FACTORY_H
#define PHASE_SOLVER_FACTORY_H

#include "Solver.h"

class SolverFactory {
public:
  enum Type { POISSON };

  static std::shared_ptr<Solver>
  create(Type type, const Input &input,
         const std::shared_ptr<StructuredGrid2D> &grid);

  static std::shared_ptr<Solver>
  create(std::string type, const Input &input,
         const std::shared_ptr<StructuredGrid2D> &grid);

  static std::shared_ptr<Solver>
  create(const Input &input, const std::shared_ptr<StructuredGrid2D> &grid);
};

#endif
