#ifndef PISO_H
#define PISO_H

#include "Solver.h"
#include "Communicator.h"
#include "FiniteVolumeEquation.h"
#include "ForceIntegrator.h"
#include "ScalarGradient.h"

class Piso : public Solver
{
public:

    Piso(const Input &input,
         std::shared_ptr<FiniteVolumeGrid2D>& grid);

    virtual Scalar solve(Scalar timeStep);

    Scalar maxCourantNumber(Scalar timeStep) const;

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    VectorFiniteVolumeField &u;
    ScalarGradient &gradP, &gradPCorr;
    ScalarFiniteVolumeField &p, &pCorr, &rho, &mu, &m, &d;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePCorrEqn();

    virtual void rhieChowInterpolation();

    void correctVelocity();

    Equation<Vector2D> uEqn_;
    Equation<Scalar> pCorrEqn_;

    size_t nInnerIterations_, nPCorrections_;
    Scalar momentumOmega_, pCorrOmega_;

    CellZone &fluid_;
};

#endif
