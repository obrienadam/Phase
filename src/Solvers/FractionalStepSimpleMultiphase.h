#ifndef FRACTIONAL_STEP_SIMPLE_MULTIPHASE_H
#define FRACTIONAL_STEP_SIMPLE_MULTIPHASE_H

#include "FractionalStepSimple.h"
#include "Celeste.h"

class FractionalStepSimpleMultiphase: public FractionalStepSimple
{
public:
    FractionalStepSimpleMultiphase(const Input& input,
                                   const Communicator& comm,
                                   std::shared_ptr<FiniteVolumeGrid2D> &grid);

    void initialize();

    Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &rho, &mu, &gamma;
    VectorFiniteVolumeField &rhoU, &gradGamma, &ft;

private:

    Scalar solveGammaEqn(Scalar timeStep);

    Scalar solveUEqn(Scalar timeStep);

    Scalar solvePEqn(Scalar timeStep);

    void correctVelocity(Scalar timeStep);

    void updateProperties(Scalar timeStep);

    //- Properties
    Scalar rho1_, rho2_, mu1_, mu2_;
    Vector2D g_;

    //- Models
    Celeste surfaceTensionModel_;

    //- Equations
    Equation<Scalar> gammaEqn_;
};

#endif