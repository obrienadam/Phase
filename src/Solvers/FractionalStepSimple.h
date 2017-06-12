#ifndef FRACTIONAL_STEP_SIMPLE_H
#define FRACTIONAL_STEP_SIMPLE_H

#include "Solver.h"
#include "FiniteVolumeEquation.h"

class FractionalStepSimple: public Solver
{
public:

    FractionalStepSimple(const Input& input,
                         const Communicator& comm,
                         std::shared_ptr<FiniteVolumeGrid2D> &grid);

    void initialize();

    std::string info() const;

    Scalar solve(Scalar timeStep);

    Scalar maxCourantNumber(Scalar timeStep) const;

    Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    VectorFiniteVolumeField &u, &gradP;
    ScalarFiniteVolumeField &p;

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
