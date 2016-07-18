#ifndef PISO_Multiphase_H
#define PISO_Multiphase_H

#include "Piso.h"
#include "SurfaceTensionForce.h"

class PisoMultiphase : public Piso
{
public:

    enum InterfaceAdvection{CICSAM, PLIC};
    enum CurvatureEvaluation{CSF, HF, CELESTE};

    PisoMultiphase(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &gamma;
    VectorFiniteVolumeField &ft;

protected:

    virtual void computeRho();
    virtual void computeMu();

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solveGammaEqn(Scalar timeStep);

    virtual void rhieChowInterpolation();

    Scalar rho1_, rho2_, mu1_, mu2_;
    std::shared_ptr<SurfaceTensionForce> surfaceTensionForce_;

    Equation<ScalarFiniteVolumeField> gammaEqn_;

    InterfaceAdvection interfaceAdvectionMethod_;
    CurvatureEvaluation curvatureEvaluationMethod_;
};

#endif
