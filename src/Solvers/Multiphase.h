#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "Piso.h"

class Multiphase : public Piso
{
public:
    Multiphase(const FiniteVolumeGrid2D& grid, const Input& input);

    Scalar solve(Scalar timeStep);

    ScalarFiniteVolumeField &gamma, &kappa;
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
};

namespace hc
{

Scalar betaFace(Scalar gammaD, Scalar gammaA, Scalar gammaU, Scalar coD);

}

namespace sc
{

}

namespace uq
{

}

namespace cicsam
{
enum Type{HC, SC, UQ};

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field, Scalar timeStep);
}

#endif
