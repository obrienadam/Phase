#ifndef FRACTIONAL_STEP_H
#define FRACTIONAL_STEP_H

#include "Solver.h"
#include "FiniteVolumeEquation.h"
#include "ScalarGradient.h"
#include "JacobianField.h"

class FractionalStep: public Solver
{
public:

    FractionalStep(const Input& input,
                         std::shared_ptr<FiniteVolumeGrid2D> &grid);

    void initialize();

    std::string info() const;

    Scalar solve(Scalar timeStep);

    Scalar maxCourantNumber(Scalar timeStep) const;

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    VectorFiniteVolumeField &u;
    ScalarFiniteVolumeField &p;
    ScalarGradient &gradP;
    JacobianField &gradU;

protected:

    Scalar solveUEqn(Scalar timeStep);

    Scalar solvePEqn(Scalar timeStep);

    void correctVelocity(Scalar timeStep);

    Scalar maxDivergenceError();

    Equation<Vector2D> uEqn_;
    Equation<Scalar> pEqn_;

    Scalar rho_, mu_;

    CellZone &fluid_;
};

#endif
