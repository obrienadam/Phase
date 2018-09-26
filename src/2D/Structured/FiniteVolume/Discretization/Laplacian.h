#ifndef PHASE_LAPLACIAN_H
#define PHASE_LAPLACIAN_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace fv
{

template<class T>
FiniteVolumeEquation<T> lap(Scalar gamma, Field<T> &phi)
{
    FiniteVolumeEquation<T> eqn(phi);

    for(const Cell& cell: phi.grid()->localCells())
    {
        for(const auto &st: cell.faceStencils())
        {
            if(st.isBoundary())
            {
                switch(phi.bctype(st.zeta()))
                {

                }
            }
            else
            {

            }
        }
    }

    return eqn;
}

}

#endif
