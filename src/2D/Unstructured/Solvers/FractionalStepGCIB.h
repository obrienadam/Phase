#ifndef PHASE_FRACTIONAL_STEP_GHOST_CELL_H
#define PHASE_FRACTIONAL_STEP_GHOST_CELL_H

#include "FiniteVolume/ImmersedBoundary/GhostCellImmersedBoundary.h"
#include "FractionalStep.h"

class FractionalStepGCIB : public FractionalStep {
public:
  FractionalStepGCIB(const Input &input,
                     const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

  Scalar solve(Scalar timeStep) override;

protected:
  Scalar solveUEqn(Scalar timeStep) override;

  Scalar solvePEqn(Scalar timeStep) override;

  GhostCellImmersedBoundary ib_;
};

#endif
