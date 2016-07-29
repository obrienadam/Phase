#ifndef FRACTIONAL_STEP_H
#define FRACTIONAL_STEP_H

#include "Solver.h"
#include "FiniteVolumeEquation.h"

class FractionalStep: public Solver
{
public:

    FractionalStep(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual std::string info() const;

    virtual Scalar solve(Scalar timeStep);
    virtual Scalar computeMaxTimeStep(Scalar maxCo) const;

    VectorFiniteVolumeField &u, &sg, &gradP, &gradSp;
    ScalarFiniteVolumeField &p, &sp, &rho, &mu, &divUStar;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solvePEqn(Scalar timeStep);
    virtual void correctVelocity(Scalar timeStep);

    void computeAdvectingVelocity(Scalar timeStep);

    Scalar courantNumber(Scalar timeStep);

    Vector2D g_;

    Equation<VectorFiniteVolumeField> uEqn_;
    Equation<ScalarFiniteVolumeField> pEqn_;

};

#endif
