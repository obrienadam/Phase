#ifndef PISO_Multiphase_H
#define PISO_Multiphase_H

#include "Piso.h"
#include "SurfaceTensionForce.h"

class PisoMultiphase : public Piso
{
public:

    enum InterfaceAdvection{CICSAM, PLIC};
    enum CurvatureEvaluation{CSF, HF, CELESTE};

    PisoMultiphase(const Input& input, const Communicator &comm, FiniteVolumeGrid2D &grid);

    virtual Scalar solve(Scalar timeStep);
    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const;

    ScalarFiniteVolumeField &gamma;
    VectorFiniteVolumeField &gradGamma, &ft, &gradRho, &sg;

protected:

    virtual void computeRho();
    virtual void computeMu();

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual void rhieChowInterpolation();

    Vector2D g_;
    Scalar rho1_, rho2_, mu1_, mu2_;
    Scalar capillaryTimeStep_;

    std::shared_ptr<SurfaceTensionForce> surfaceTensionForce_;

    Equation<ScalarFiniteVolumeField> gammaEqn_;

    InterfaceAdvection interfaceAdvectionMethod_;
    CurvatureEvaluation curvatureEvaluationMethod_;
};

#endif
