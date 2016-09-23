#ifndef COUPLED_H
#define COUPLED_H

#include "Solver.h"
#include "CoupledEquation.h"

class Coupled: public Solver
{
public:

    Coupled(const Input& input, FiniteVolumeGrid2D &grid);

    virtual Scalar solve(Scalar timeStep);
    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    VectorFiniteVolumeField &u;
    ScalarFiniteVolumeField &p, &rho, &mu;
    CoupledEquation eqn_;

protected:

    Vector2D g_;

};

#endif
