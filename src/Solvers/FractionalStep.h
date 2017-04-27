#ifndef FRACTIONAL_STEP_H
#define FRACTIONAL_STEP_H

#include <fstream>

#include "Solver.h"
#include "Communicator.h"
#include "FiniteVolumeEquation.h"
#include "ForceIntegrator.h"
#include "TecplotViewer.h"

class FractionalStep : public Solver
{
public:

    FractionalStep(const Input &input, const Communicator &comm, FiniteVolumeGrid2D &grid);
    ~FractionalStep() {log.close();}

    virtual void initialize();

    virtual std::string info() const;

    virtual Scalar solve(Scalar timeStep);

    Scalar maxCourantNumber(Scalar timeStep) const;

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    VectorFiniteVolumeField &u, &gradP;
    ScalarFiniteVolumeField &p, &rho, &mu;

protected:

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solvePEqn(Scalar timeStep);

    virtual void correctVelocity(Scalar timeStep);

    virtual void computeFaceVelocities(Scalar timeStep);

    Vector2D maxVelocity() const;

    Vector2D maxFaceVelocity() const;

    Scalar maxDivergenceError() const;

    Scalar alphaAdv_, alphaDiff_;

    Equation<Vector2D> uEqn_;
    Equation<Scalar> pEqn_;

    std::ofstream log;
};

#endif
