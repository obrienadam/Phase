#ifndef PHASE_FRACTIONAL_STEP_H
#define PHASE_FRACTIONAL_STEP_H

#include "FiniteVolume/Field/ScalarGradient.h"
#include "FiniteVolume/Field/JacobianField.h"

#include "Solver.h"

class FractionalStep: public Solver
{
public:

    FractionalStep(const Input& input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    virtual void initialize();

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

    FiniteVolumeEquation<Vector2D> uEqn_;

    FiniteVolumeEquation<Scalar> pEqn_;

    Scalar rho_, mu_;

    Vector2D g_;

    CellGroup &fluid_;
};

#endif
