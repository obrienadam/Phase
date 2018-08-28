#ifndef PHASE_LAPLACIAN_H
#define PHASE_LAPLACIAN_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"
#include "Math/TaylorSeries.h"

namespace fv
{

template<class T>
FiniteVolumeEquation<T> lap(Scalar gamma, Field<T> &phi)
{
    FiniteVolumeEquation<T> eqn(phi);

    for(const Cell& cell: phi.grid()->cells())
    {
        for(auto dir: Stencil::DIRECTIONS)
        {
            if(cell.face(dir).isBoundary())
            auto st = Stencil(cell, {1, 0}, dir);

        }
    }

    return eqn;
}

}

#endif
