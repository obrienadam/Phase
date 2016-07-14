#ifndef FRACTIONAL_STEP_H
#define FRACTIONAL_STEP_H

#include "Solver.h"
#include "FiniteVolumeEquation.h"

class FractionalStep: public Solver
{
public:

    FractionalStep(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual Scalar solve(Scalar timeStep);
    virtual Scalar computeMaxTimeStep(Scalar maxCo) const;

    VectorFiniteVolumeField &u;
    ScalarFiniteVolumeField &p, &rho, &mu, &divUStar;

protected:

    Scalar solveUEqn(Scalar timeStep);
    Scalar solvePEqn(Scalar timeStep);
    void correctVelocity(Scalar timeStep);

    Equation<VectorFiniteVolumeField> uEqn_;
    Equation<ScalarFiniteVolumeField> pEqn_;
};

#endif
