#ifndef PHASE_SOLVER_FACTORY_H
#define PHASE_SOLVER_FACTORY_H

#include "Solver.h"

class SolverFactory {
public:
  enum SolverType { POISSON, FRACTIONAL_STEP };

  static std::shared_ptr<Solver>
  create(SolverType type, const Input &input,
         const std::shared_ptr<const StructuredGrid3D> &grid);

  static std::shared_ptr<Solver>
  create(std::string type, const Input &input,
         const std::shared_ptr<const StructuredGrid3D> &grid);

  static std::shared_ptr<Solver>
  create(const Input &input,
         const std::shared_ptr<const StructuredGrid3D> &grid);

private:
};

#endif
