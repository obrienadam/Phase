#ifndef PHASE_FRACTIONAL_STEP_BOUSSINESQ_H
#define PHASE_FRACTIONAL_STEP_BOUSSINESQ_H

#include "FractionalStep.h"

class FractionalStepBoussinesq : public FractionalStep {
public:
  FractionalStepBoussinesq(
      const Input &input,
      const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

  Scalar solve(Scalar timeStep);

  ScalarFiniteVolumeField &T;

protected:
  FiniteVolumeEquation<Scalar> TEqn_;

  Scalar solveUEqn(Scalar timeStep);

  Scalar solveTEqn(Scalar timeStep);

  Scalar alpha_, T0_, kappa_;
};

#endif
