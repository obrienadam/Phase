#ifndef COUPLED_H
#define COUPLED_H

#include "Solver.h"
#include "CoupledEquation.h"

class Coupled: public Solver
{
public:

    Coupled(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual Scalar solve(Scalar timeStep);
    virtual Scalar computeMaxTimeStep(Scalar maxCo) const;

    VectorFiniteVolumeField &u;
    ScalarFiniteVolumeField &p, &rho, &mu;
    CoupledEquation eqn_;

protected:

    Vector2D g_;

};

#endif
