#ifndef FRACTIONAL_STEP_H
#define FRACTIONAL_STEP_H

#include "Solver.h"
#include "ScalarGradient.h"
#include "JacobianField.h"

class FractionalStep: public Solver
{
public:

    FractionalStep(const Input& input);

    void initialize();

    std::string info() const;

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

    virtual Scalar maxDivergenceError();

    Equation<Vector2D> uEqn_;
    Equation<Scalar> pEqn_;

    Scalar rho_, mu_;
    Vector2D g_;

    CellZone &fluid_;
};

#endif
