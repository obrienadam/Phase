#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "Piso.h"
#include "SurfaceTensionForce.h"

class Multiphase : public Piso
{
public:

    enum InterfaceAdvection{CICSAM, PLIC};
    enum CurvatureEvaluation{CSF, HF};

    Multiphase(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &gamma;
    VectorFiniteVolumeField &ft;

protected:

    virtual void computeRho();
    virtual void computeMu();

    virtual Scalar solveUEqn(Scalar timeStep);
    virtual Scalar solveGammaEqn(Scalar timeStep);

    Scalar rho1_, rho2_, mu1_, mu2_;
    std::shared_ptr<SurfaceTensionForce> surfaceTensionForce_;

    Equation<ScalarFiniteVolumeField> gammaEqn_;

    InterfaceAdvection interfaceAdvectionMethod_;
    CurvatureEvaluation curvatureEvaluationMethod_;
};

#endif
