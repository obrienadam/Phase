#ifndef PISO_H
#define PISO_H

#include "Solver.h"
#include "Input.h"
#include "FiniteVolumeEquation.h"

class Piso : public Solver
{
public:

    Piso(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual Scalar solve(Scalar timeStep);
    Scalar computeMaxTimeStep(Scalar maxCo) const;

    VectorFiniteVolumeField &u, &h, &sg, &gradP, &gradPCorr;
    ScalarFiniteVolumeField &p, &pCorr, &rho, &mu, &m, &d;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solvePCorrEqn();

    virtual void rhieChowInterpolation();

    void correctAll();
    void correctPressure();
    void correctVelocity();

    void computeStaticPressure();

    Scalar courantNumber(Scalar timeStep);

    Vector2D g_;

    Equation<VectorFiniteVolumeField> uEqn_;
    Equation<ScalarFiniteVolumeField> pCorrEqn_;

    size_t nInnerIterations_, nPCorrections_;
    Scalar momentumOmega_, pCorrOmega_;
};

#endif
