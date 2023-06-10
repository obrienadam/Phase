#ifndef PHASE_POISSON_H
#define PHASE_POISSON_H

#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"

#include "Solver.h"

class Poisson : public Solver {
public:
  Poisson(const Input &input,
          const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

  void initialize();

  Scalar solve(Scalar timeStep);

  Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const {
    return std::numeric_limits<Scalar>::infinity();
  }

  ScalarFiniteVolumeField &phi;

protected:
  FiniteVolumeEquation<Scalar> phiEqn_;

  std::shared_ptr<CellGroup> solid_;

  Scalar gamma_;
};

#endif
