#ifndef CICSAM_H
#define CICSAM_H

#include "Equation.h"

namespace cicsam
{
    ScalarFiniteVolumeField beta(const VectorFiniteVolumeField &u,
                                 const VectorFiniteVolumeField &gradGamma,
                                 const ScalarFiniteVolumeField &gamma,
                                 Scalar timeStep,
                                 Scalar k = 1.);

    Equation<Scalar> div(const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &beta,
                         ScalarFiniteVolumeField &gamma,
                         const CellGroup &cells,
                         Scalar theta);

    Equation<Scalar> div(const VectorFiniteVolumeField &u,
                         const ScalarFiniteVolumeField &beta,
                         ScalarFiniteVolumeField &gamma,
                         Scalar theta);
}

#endif
