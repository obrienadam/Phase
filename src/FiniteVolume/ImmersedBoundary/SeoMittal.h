#ifndef SEO_MITTAL_H
#define SEO_MITTAL_H

#include "ImmersedBoundary.h"
#include "FiniteVolumeEquation.h"

namespace seo
{

    //- Constant density
    Equation<Scalar> laplacian(const ImmersedBoundary &ib, Scalar rho, Scalar timeStep, ScalarFiniteVolumeField& p);

    ScalarFiniteVolumeField div(const ImmersedBoundary &ib, const VectorFiniteVolumeField& u);

    void correct(const ImmersedBoundary &ib, Scalar rho, const ScalarFiniteVolumeField& p,
                 VectorFiniteVolumeField& gradP, VectorFiniteVolumeField& u, Scalar timeStep);

    //- Variable density, e.g. multiphase
    Equation<Scalar> laplacian(const ImmersedBoundary &ib, const ScalarFiniteVolumeField& rho, Scalar timeStep,
                               ScalarFiniteVolumeField &p);

    void correct(const ImmersedBoundary &ib, const ScalarFiniteVolumeField &rho, const ScalarFiniteVolumeField& p,
                 VectorFiniteVolumeField &gradP, VectorFiniteVolumeField &u, Scalar timeStep);

    //- Divergence
    Scalar maxDivergence(const ImmersedBoundary &ib, const VectorFiniteVolumeField& u);
};

#endif