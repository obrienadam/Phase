#ifndef PHASE_LAPLACIAN_H
#define PHASE_LAPLACIAN_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace fv
{

template<class T>
FiniteVolumeEquation<T> lap(Scalar gamma, Field<T> &phi, Scalar theta = 1.)
{
    FiniteVolumeEquation<T> eqn(phi);
    const Field<T> &phi0 = phi.oldField(0);

    for(const Cell& cell: phi.grid()->localCells())
    {
        for(const InteriorFaceStencil &st: cell.interiorFaceStencils())
        {
            Scalar an = gamma * dot(st.rc(1), st.sf()) / st.rc(1).magSqr();

            eqn.add(cell, cell, -an * theta);
            eqn.add(cell, st.cell(1), an * theta);
            eqn.addSource(cell, an * (phi0(st.cell(1)) - phi0(st.cell(0))) * (1. - theta));
        }
    }

    return eqn;
}

}

#endif
