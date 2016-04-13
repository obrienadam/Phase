#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "Piso.h"

class Multiphase : public Piso
{
public:
    Multiphase(const FiniteVolumeGrid2D& grid, const Input& input);

    Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &gamma;

private:

    void computeRho();
    void computeMu();
    Scalar solveGammaEqn(Scalar timeStep);

    Scalar rho1_, rho2_, mu1_, mu2_, sigma_;

    Equation<ScalarFiniteVolumeField> gammaEqn_;
};

#endif
