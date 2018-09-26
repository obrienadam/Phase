#ifndef PHASE_AXISYMMETRIC_TIME_DERIVATIVE_H
#define PHASE_AXISYMMETRIC_TIME_DERIVATIVE_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace axi
{
    template<class T>
    FiniteVolumeEquation<T> ddt(FiniteVolumeField<T> &phi, Scalar timeStep)
    {
        FiniteVolumeEquation<T> eqn(phi);
        const FiniteVolumeField<T> &phi0 = phi.oldField(0);

        for (const Cell &cell: phi.cells())
        {
            Scalar volume = cell.polarVolume();
            eqn.add(cell, cell, volume / timeStep);
            eqn.addSource(cell, -phi0(cell) * volume / timeStep);
        }

        return eqn;
    }
}

#endif
