#ifndef AXISYMMETRIC_TIME_DERIVATIVE_H
#define AXISYMMETRIC_TIME_DERIVATIVE_H

#include "Equation.h"

namespace axi
{
    template<class T>
    Equation<T> ddt(FiniteVolumeField<T> &phi,
                    Scalar timeStep)
    {
        Equation<T> eqn(phi);

        auto polarVolume = [](const Cell &cell) {
            Scalar volume = 0.;

            for (const InteriorLink &nb: cell.neighbours())
                volume += dot(nb.face().centroid(), nb.face().polarOutwardNorm(cell.centroid()));

            for (const BoundaryLink &bd: cell.boundaries())
                volume += dot(bd.face().centroid(), bd.face().polarOutwardNorm(cell.centroid()));

            return volume / 3.;
        };

        for (const Cell &cell: phi.cells())
        {
            Scalar volume = polarVolume(cell);
            eqn.add(cell, cell, volume / timeStep);
            eqn.addSource(cell, -phi(cell) * volume / timeStep);
        }

        return eqn;
    }
}

#endif
