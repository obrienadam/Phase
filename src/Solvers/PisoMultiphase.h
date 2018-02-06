#ifndef PISO_MULTIPHASE_H
#define PISO_MULTIPHASE_H

#include "Piso.h"
#include "Multiphase/SurfaceTensionForce.h"

class PisoMultiphase : public Piso
{
public:

    enum InterfaceAdvection
    {
        CICSAM
    };

    PisoMultiphase(const Input &input,
                   std::shared_ptr<FiniteVolumeGrid2D>& grid);

    virtual Scalar solve(Scalar timeStep);

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    ScalarFiniteVolumeField &gamma, &beta;
    ScalarGradient &gradGamma, &gradRho;
    VectorFiniteVolumeField &sg;

protected:

    virtual void computeRho();

    virtual void computeMu();

    virtual Scalar solveUEqn(Scalar timeStep);

    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual void rhieChowInterpolation();

    Vector2D g_;
    Scalar rho1_, rho2_, mu1_, mu2_;
    Scalar capillaryTimeStep_;

    std::shared_ptr<SurfaceTensionForce> ft_;

    Equation<Scalar> gammaEqn_;

    InterfaceAdvection interfaceAdvectionMethod_;
    //CurvatureEvaluation curvatureEvaluationMethod_;
};

#endif
