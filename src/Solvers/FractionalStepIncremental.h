#ifndef FRACTIONAL_STEP_INCREMENTAL_H
#define FRACTIONAL_STEP_INCREMENTAL_H

#include <fstream>

#include "Solver.h"
#include "FiniteVolumeEquation.h"
#include "ScalarGradient.h"
#include "JacobianField.h"

class FractionalStepIncremental : public Solver
{
public:

    FractionalStepIncremental(const Input &input,
                   std::shared_ptr<FiniteVolumeGrid2D> &grid);

    virtual void initialize();

    virtual std::string info() const;

    virtual Scalar solve(Scalar timeStep);

    Scalar maxCourantNumber(Scalar timeStep) const;

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    VectorFiniteVolumeField &u;
    ScalarFiniteVolumeField &p;
    ScalarGradient &gradP;
    JacobianField &gradU;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    virtual void correctVelocity(Scalar timeStep);

    Scalar maxDivergenceError() const;

    Scalar rho_, mu_;

    Equation<Vector2D> uEqn_;
    Equation<Scalar> pEqn_;

    CellZone& fluid_;
};

#endif
