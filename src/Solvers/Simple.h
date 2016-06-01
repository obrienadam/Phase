#ifndef SIMPLE_H
#define SIMPLE_H

#include "Solver.h"
#include "Input.h"
#include "FiniteVolumeEquation.h"

class Simple : public Solver
{
public:

    Simple(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual Scalar solve(Scalar timeStep);
    Scalar computeMaxTimeStep(Scalar maxCo) const;

    VectorFiniteVolumeField &u, &h;
    ScalarFiniteVolumeField &p, &pCorr, &rho, &mu, &m, &d, &ids, &globalIndices;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solvePCorrEqn();

    void computeD();
    void computeH();
    void rhieChowInterpolation();

    void correctAll();
    void correctPressure();
    void correctVelocity();

    Scalar courantNumber(Scalar timeStep);

    Vector2D g_;

    Equation<VectorFiniteVolumeField> uEqn_;
    Equation<ScalarFiniteVolumeField> pCorrEqn_;

    size_t nInnerIterations_;
    Scalar momentumOmega_, pCorrOmega_;
};

#endif
