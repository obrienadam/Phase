#ifndef FRACTIONAL_STEP_H
#define FRACTIONAL_STEP_H

#include "Solver.h"
#include "Communicator.h"
#include "FiniteVolumeEquation.h"
#include "ForceIntegrator.h"

class FractionalStep: public Solver
{
public:

    FractionalStep(const Input& input, const Communicator& comm, FiniteVolumeGrid2D &grid);

    virtual std::string info() const;

    virtual Scalar solve(Scalar timeStep);
    Scalar maxCourantNumber(Scalar timeStep) const;
    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    VectorFiniteVolumeField &u, &gradP, &gradDp;
    ScalarFiniteVolumeField &p, &dp, &rho, &mu, &divUStar;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solvePEqn(Scalar timeStep);
    virtual void correctVelocity(Scalar timeStep);

    virtual void computeFaceVelocities(Scalar timeStep);
    virtual void computeMassSource(Scalar timeStep);

    Equation<Vector2D> uEqn_;
    Equation<Scalar> pEqn_;
};

#endif
