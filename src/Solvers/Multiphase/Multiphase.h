#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "Piso.h"
#include "RangeSearch.h"

class Multiphase : public Piso
{
public:
    Multiphase(const FiniteVolumeGrid2D& grid, const Input& input);

    Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &gamma, &gammaTilde, &kappa;
    VectorFiniteVolumeField &n;

protected:

    void computeRho();
    void computeMu();

    Scalar solveUEqn(Scalar timeStep);
    Scalar solveGammaEqn(Scalar timeStep);

    void computeInterfaceNormals();
    void computeCurvature();

    Scalar rho1_, rho2_, mu1_, mu2_, sigma_;

    Equation<ScalarFiniteVolumeField> gammaEqn_;

    RangeSearch cellRangeSearch_;
    Scalar kernelWidth_;
};

#endif
