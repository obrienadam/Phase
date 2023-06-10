#ifndef PHASE_COUPLED_FINITE_VOLUME_EQUATION_H
#define PHASE_COUPLED_FINITE_VOLUME_EQUATION_H

#include "FiniteVolumeEquation.h"
#include "Math/Equation.h"

class CoupledEquation : public Equation {
public:
  void add(const FiniteVolumeEquation<Vector2D> &uEqn,
           const FiniteVolumeEquation<Scalar> &pEqn);

  void add(const FiniteVolumeEquation<Scalar> &pEqn,
           const FiniteVolumeEquation<Vector2D> &uEqn);

protected:
};

#endif
