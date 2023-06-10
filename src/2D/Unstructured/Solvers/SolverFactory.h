#ifndef PHASE_SOLVER_FACTORY_H
#define PHASE_SOLVER_FACTORY_H

#include "Solver.h"

class SolverFactory {
public:
  enum SolverType {
    POISSON,
    FRACTIONAL_STEP,
    FRACTIONAL_STEP_AXISYMMETRIC,
    FRACTIONAL_STEP_AXISYMMETRIC_DIRECT_FORCING_IB,
    FRACTIONAL_STEP_AXISYMMETRIC_DIRECT_FORCING_IB_MULTIPHASE,
    FRACTIONAL_STEP_MULTIPHASE,
    FRACTIONAL_STEP_EULER_LAGRANGE_IB,
    FRACTIONAL_STEP_DIRECT_FORCING_IB,
    FRACTIONAL_STEP_DIRECT_FORCING_IB_MULTIPHASE,
    FRACTIONAL_STEP_GHOST_CELL_IB,
    FRACTIONAL_STEP_GHOST_CELL_IB_MULTIPHASE,
    FRACTIONAL_STEP_BOUSSINESQ
  };

  static std::shared_ptr<Solver>
  create(SolverType type, const Input &input,
         const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

  static std::shared_ptr<Solver>
  create(std::string type, const Input &input,
         const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

  static std::shared_ptr<Solver>
  create(const Input &input,
         const std::shared_ptr<const FiniteVolumeGrid2D> &grid);
};

#endif // PHASE_SOLVERFACTORY_H
