#ifndef PHASE_LAPLACIAN_H
#define PHASE_LAPLACIAN_H

#include "FiniteVolumeEquation.h"

namespace fv
{
    template<class T>
    FiniteVolumeEquation<T> laplacian(Scalar gamma, FiniteVolumeField<T> &phi, Scalar theta = 1.)
    {
        FiniteVolumeEquation<T> eqn(phi);

        for (const Cell &cell: phi.cells())
        {
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar coeff = gamma * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
                eqn.add(cell, nb.cell(), theta * coeff);
                eqn.add(cell, cell, theta * -coeff);
                eqn.addSource(cell, (1. - theta) * coeff * (phi(nb.cell()) - phi(cell)));
            }

            for (const BoundaryLink &bd: cell.boundaries())
            {
                Scalar coeff = gamma * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();

                switch (phi.boundaryType(bd.face()))
                {
                    case FiniteVolumeField<T>::FIXED:
                        eqn.add(cell, cell, theta * -coeff);
                        eqn.addSource(cell, theta * coeff * phi(bd.face()));
                        eqn.addSource(cell, (1. - theta) * coeff * (phi(bd.face()) - phi(cell)));
                        break;

                    case FiniteVolumeField<T>::NORMAL_GRADIENT:
                    case FiniteVolumeField<T>::SYMMETRY:
                        break;

                    default:
                        throw Exception("fv", "laplacian<T>", "unrecognized or unspecified boundary type.");
                }
            }
        }

        return eqn;
    }

    template<class T>
    FiniteVolumeEquation<T> laplacian(const ScalarFiniteVolumeField &gamma,
                          FiniteVolumeField<T> &phi,
                          Scalar theta = 1.)
    {
        FiniteVolumeEquation<T> eqn(phi);
        const ScalarFiniteVolumeField &gamma0 = gamma.oldField(0);

        for (const Cell &cell: phi.cells())
        {
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar coeff = gamma(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
                Scalar coeff0 = gamma0(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) / nb.rCellVec().magSqr();
                eqn.add(cell, cell, theta * -coeff);
                eqn.add(cell, nb.cell(), theta * coeff);
                eqn.addSource(cell, (1. - theta) * coeff0 * (phi(nb.cell()) - phi(cell)));
            }

            for (const BoundaryLink &bd: cell.boundaries())
            {
                Scalar coeff = gamma(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();
                Scalar coeff0 = gamma0(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) / bd.rFaceVec().magSqr();

                switch (phi.boundaryType(bd.face()))
                {
                    case FiniteVolumeField<T>::FIXED:
                        eqn.add(cell, cell, theta * -coeff);
                        eqn.addSource(cell, theta * coeff * phi(bd.face()));
                        eqn.addSource(cell, (1. - theta) * coeff0 * (phi(bd.face()) - phi(cell)));
                        break;

                    case FiniteVolumeField<T>::NORMAL_GRADIENT:
                    case FiniteVolumeField<T>::SYMMETRY:
                        break;

                    default:
                        throw Exception("fv", "laplacian<T>", "unrecognized or unspecified boundary type.");
                }
            }
        }

        return eqn;
    }
}

#endif
