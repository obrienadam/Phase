#ifndef PHASE_CICSAM_H
#define PHASE_CICSAM_H

#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace cicsam
{
Scalar hc(Scalar gammaDTilde, Scalar coD);

Scalar uq(Scalar gammaDTilde, Scalar coD);

Scalar beta(const VectorFiniteVolumeField &u,
            const ScalarFiniteVolumeField &gamma,
            const VectorFiniteVolumeField &gradGamma,
            Scalar timeStep,
            Scalar k,
            const Face &f);

void computeMomentumFlux(Scalar rho1,
                         Scalar rho2,
                         const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &gamma,
                         const VectorFiniteVolumeField &gradGamma,
                         Scalar timeStep,
                         VectorFiniteVolumeField &rhoU);

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 ScalarFiniteVolumeField &gamma,
                                 const VectorFiniteVolumeField &gradGamma,
                                 Scalar timeStep,
                                 Scalar theta,
                                 const CellGroup &cells);

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 ScalarFiniteVolumeField &gamma,
                                 const VectorFiniteVolumeField &gradGamma,
                                 Scalar timeStep,
                                 Scalar theta);

FiniteVolumeEquation<Scalar> div2e(const VectorFiniteVolumeField &u,
                                   ScalarFiniteVolumeField &gamma,
                                   const VectorFiniteVolumeField &gradGamma,
                                   Scalar timeStep,
                                   Scalar theta,
                                   const CellGroup &cells);

FiniteVolumeEquation<Scalar> div2e(const VectorFiniteVolumeField &u,
                                   ScalarFiniteVolumeField &gamma,
                                   const VectorFiniteVolumeField &gradGamma,
                                   Scalar timeStep,
                                   Scalar theta);
}

#endif
