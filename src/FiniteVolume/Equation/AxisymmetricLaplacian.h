#ifndef AXISYMMETRIC_LAPLACIAN_H
#define AXISYMMETRIC_LAPLACIAN_H

#include "Equation.h"

namespace axi
{
    template<typename T>
    Equation<T> laplacian(Scalar gamma, FiniteVolumeField<T> &phi, Scalar theta = 1.)
    {
        Equation<T> eqn(phi);

        for (const Cell &cell: phi.cells())
        {
            for (const InteriorLink &nb: cell.neighbours())
            {
                Vector2D sf = nb.face().polarOutwardNorm(cell.centroid());
                Scalar coeff = theta * gamma * dot(nb.rCellVec(), sf) / nb.rCellVec().magSqr();
                Scalar coeff0 = (1. - theta) * gamma * dot(nb.rCellVec(), sf) / nb.rCellVec().magSqr();

                eqn.add(cell, cell, -coeff);
                eqn.add(cell, nb.cell(), coeff);
                eqn.addSource(cell, coeff0 * (phi(nb.cell()) - phi(cell)));
            }

            for (const BoundaryLink &bd: cell.boundaries())
            {
                Vector2D sf = bd.face().polarOutwardNorm(cell.centroid());
                Scalar coeff = theta * gamma * dot(bd.rFaceVec(), sf) / bd.rFaceVec().magSqr();
                Scalar coeff0 = (1. - theta) * gamma * dot(bd.rFaceVec(), sf) / bd.rFaceVec().magSqr();

                switch (phi.boundaryType(bd.face()))
                {
                    case FiniteVolumeField<T>::FIXED:
                        eqn.add(cell, cell, -coeff);
                        eqn.addSource(cell, coeff * phi(bd.face()));
                        eqn.addSource(cell, coeff0 * (phi(bd.face()) - phi(cell)));
                        break;

                    case FiniteVolumeField<T>::NORMAL_GRADIENT:
                    case FiniteVolumeField<T>::SYMMETRY:
                        break;

                    default:
                        throw Exception("axi", "laplacian<T>", "unrecognized or unspecified boundary type.");
                }
            }
        }

        return eqn;
    }
}

#endif
