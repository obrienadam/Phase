#ifndef AXISYMMETRIC_TIME_DERIVATIVE_H
#define AXISYMMETRIC_TIME_DERIVATIVE_H

#include "Equation.h"

namespace axi
{
    template<class T>
    Equation<T> ddt(FiniteVolumeField<T> &phi, Scalar timeStep)
    {
        Equation<T> eqn(phi);

        for (const Cell &cell: phi.cells())
        {
            Scalar volume = cell.polarVolume();
            eqn.add(cell, cell, volume / timeStep);
            eqn.addSource(cell, -phi(cell) * volume / timeStep);
        }

        return eqn;
    }
}

#endif
