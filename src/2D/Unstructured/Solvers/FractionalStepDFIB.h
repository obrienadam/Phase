#ifndef PHASE_FRACTIONAL_STEP_DIRECT_FORCING_H
#define PHASE_FRACTIONAL_STEP_DIRECT_FORCING_H

#include "FractionalStep.h"

class DirectForcingImmersedBoundary;

class FractionalStepDFIB : public FractionalStep {
public:
  FractionalStepDFIB(const Input &input,
                     const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

  virtual void initialize() override;

  virtual Scalar solve(Scalar timeStep) override;

  virtual std::shared_ptr<const ImmersedBoundary> ib() const override;

protected:
  virtual Scalar solveUEqn(Scalar timeStep) override;

  virtual void solveExtEqns();

  void computIbForce(Scalar timeStep);

  VectorFiniteVolumeField &fb_;

  FiniteVolumeEquation<Vector2D> fbEqn_, extEqn_;

  std::shared_ptr<DirectForcingImmersedBoundary> ib_;
};

#endif
