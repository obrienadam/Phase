#include "Laplacian.h"

namespace fv
{

template<class T>
FiniteVolumeEquation<T> lap(Scalar gamma, Field<T> &phi)
{
    FiniteVolumeEquation<T> eqn(phi);

    for(const Cell &cell: phi.grid()->cells())
    {
        for(const FaceStencil &nb: cell.interiorStencils())
        {
            Scalar faceArea = nb.face().norm().mag();

            for(int i = 0; i < nb.cells().size(); ++i)
                eqn.add(cell, nb.cells()[i], nb.coeffs()[i] * gamma * faceArea);
        }

        for(const BoundaryFaceStencil &bd: cell.boundaryStencils())
        {
            Scalar faceArea = bd.face().norm().mag();
        }
    }

    return eqn;
}

}
