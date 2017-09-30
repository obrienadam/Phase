#include "SeoMittal.h"
#include "GhostCellImmersedBoundaryObject.h"

Equation<Scalar> seo::laplacian(const ImmersedBoundary &ib, Scalar rho, Scalar timeStep, ScalarFiniteVolumeField& p)
{
    Equation<Scalar> pEqn(p, "pEqn");

    return pEqn;
}

void seo::correct(const ImmersedBoundary &ib, Scalar rho, const ScalarFiniteVolumeField &p, VectorFiniteVolumeField &gradP,
                  VectorFiniteVolumeField &u, Scalar timeStep)
{

}

Equation<Scalar> seo::laplacian(const ImmersedBoundary &ib, const ScalarFiniteVolumeField& rho, Scalar timeStep,
                                ScalarFiniteVolumeField &p)
{

}


void seo::correct(const ImmersedBoundary &ib,
                  const ScalarFiniteVolumeField &rho,
                  const ScalarFiniteVolumeField& p,
                  const VectorFiniteVolumeField& src,
                  VectorFiniteVolumeField &gradP,
                  VectorFiniteVolumeField &u, Scalar timeStep)
{

}

ScalarFiniteVolumeField seo::div(const ImmersedBoundary &ib, const VectorFiniteVolumeField &u)
{

}

Scalar seo::maxDivergence(const ImmersedBoundary &ib, const VectorFiniteVolumeField& u)
{
    Scalar maxDivU = 0.;

    return maxDivU;
}