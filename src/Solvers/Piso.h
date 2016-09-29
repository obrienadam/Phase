#ifndef PISO_H
#define PISO_H

#include "Solver.h"
#include "Input.h"
#include "FiniteVolumeEquation.h"
#include "ForceIntegrator.h"

class Piso : public Solver
{
public:

    Piso(const Input& input, FiniteVolumeGrid2D &grid);

    virtual Scalar solve(Scalar timeStep);

    Scalar maxCourantNumber(Scalar timeStep) const;
    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    VectorFiniteVolumeField &u, &gradP, &gradPCorr;
    ScalarFiniteVolumeField &p, &pCorr, &rho, &mu, &m, &d;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solvePCorrEqn();

    virtual void rhieChowInterpolation();

    void correctAll();
    void correctPressure();
    void correctVelocity();

    Equation<VectorFiniteVolumeField> uEqn_;
    Equation<ScalarFiniteVolumeField> pCorrEqn_;

    size_t nInnerIterations_, nPCorrections_;
    Scalar momentumOmega_, pCorrOmega_;
};

#endif
