#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "Piso.h"

class Multiphase : public Piso
{
public:

    enum InterfaceAdvection{CICSAM, PLIC};
    enum CurvatureEvaluation{CSF, HF};

    Multiphase(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual Scalar solve(Scalar timeStep, Scalar prevTimeStep);

    ScalarFiniteVolumeField &gamma, &gammaTilde, &kappa;
    VectorFiniteVolumeField &n;

protected:

    void constructSmoothingKernels();

    virtual void computeRho();
    virtual void computeMu();

    virtual Scalar solveUEqn(Scalar timeStep, Scalar prevTimeStep);
    virtual Scalar solveGammaEqn(Scalar timeStep, Scalar prevTimeStep);

    void computeInterfaceNormals();
    void computeCurvature();

    Scalar rho1_, rho2_, mu1_, mu2_, sigma_;

    Equation<ScalarFiniteVolumeField> gammaEqn_;

    std::vector< std::vector< Ref<const Cell> > > cellRangeSearch_;
    Scalar kernelWidth_;

    InterfaceAdvection interfaceAdvectionMethod_;
    CurvatureEvaluation curvatureEvaluationMethod_;
};

#endif
