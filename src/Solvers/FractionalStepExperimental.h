#ifndef FRACTIONAL_STEP_EXPERIMENTAL_H
#define FRACTIONAL_STEP_EXPERIMENTAL_H

#include "Solver.h"

class FractionalStepExperimental: public Solver
{
public:
    FractionalStepExperimental(const Input& input, const Communicator& comm, FiniteVolumeGrid2D &grid);

    virtual std::string info() const;

    virtual Scalar solve(Scalar timeStep);
    Scalar maxCourantNumber(Scalar timeStep) const;
    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    VectorFiniteVolumeField &u, &gradP, &gradPhi;
    ScalarFiniteVolumeField &p, &phi, &rho, &mu, &divUStar;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solvePhiEqn(Scalar timeStep);
    virtual void correctVelocity(Scalar timeStep);

    virtual void computeFaceVelocities(Scalar timeStep);
    virtual void computeMassSource(Scalar timeStep);

    Equation<Vector2D> uEqn_;
    Equation<Scalar> phiEqn_;
};

#endif
