#ifndef PHASE_FRACTIONAL_STEP_H
#define PHASE_FRACTIONAL_STEP_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"
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

    virtual Scalar maxCourantNumber(Scalar timeStep) const;

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    virtual void correctVelocity(Scalar timeStep);

    virtual Scalar maxDivergenceError();

    Scalar rho_, mu_;

    Vector2D g_;

    std::shared_ptr<CellGroup> fluid_;

    VectorFiniteVolumeField &u_;

    ScalarFiniteVolumeField &p_, &co_;

    ScalarGradient &gradP_;

    JacobianField &gradU_;

    FiniteVolumeEquation<Vector2D> uEqn_;

    FiniteVolumeEquation<Scalar> pEqn_;

};

#endif
