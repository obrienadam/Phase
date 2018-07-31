#include "Laplacian.h"

namespace fv
{

template<class T>
FiniteVolumeEquation<T> lap(Field<T> &phi)
{
    FiniteVolumeEquation<T> eqn(phi);

    for(const Cell &cell: phi.grid()->cells())
    {
        for(const FaceStencil &nb: cell.interiorStencils())
        {

        }

        for(const FaceStencil &bd: cell.boundaryStencils())
        {

        }
    }

    return eqn;
}

}
