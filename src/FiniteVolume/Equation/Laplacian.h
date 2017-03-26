#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "Equation.h"

namespace fv
{

template<typename T>
Equation<T> laplacian(const ScalarFiniteVolumeField &gamma, FiniteVolumeField<T> &field)
{
    Equation<T> eqn(field);

    for (const Cell &cell: field.grid.cellZone("fluid"))
    {
        Scalar centralCoeff = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar coeff = gamma(nb.face()) * dot(nb.rCellVec(), nb.outwardNorm()) / dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;
            eqn.add(cell, nb.cell(), coeff);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar coeff =
                    gamma(bd.face()) * dot(bd.rFaceVec(), bd.outwardNorm()) / dot(bd.rFaceVec(), bd.rFaceVec());

            switch (field.boundaryType(bd.face()))
            {
            case FiniteVolumeField<T>::FIXED:
                centralCoeff -= coeff;
                eqn.addBoundary(cell, -coeff * field(bd.face()));
                break;

            case FiniteVolumeField<T>::NORMAL_GRADIENT:
            case FiniteVolumeField<T>::SYMMETRY:
                break;

            default:
                throw Exception("fv", "laplacian<T>", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(cell, cell, centralCoeff);
    }

    return eqn;
}

template<typename T>
Equation<T> laplacian(Scalar gamma, FiniteVolumeField<T> &field)
{
    Equation<T> eqn(field);

    for (const Cell &cell: field.grid.cellZone("fluid"))
    {
        Scalar centralCoeff = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar coeff = gamma * dot(nb.rCellVec(), nb.outwardNorm()) / dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;
            eqn.add(cell, nb.cell(), coeff);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar coeff = gamma * dot(bd.rFaceVec(), bd.outwardNorm()) / dot(bd.rFaceVec(), bd.rFaceVec());

            switch (field.boundaryType(bd.face()))
            {
            case FiniteVolumeField<T>::FIXED:
                centralCoeff -= coeff;
                eqn.addBoundary(cell, -coeff * field(bd.face()));
                break;

            case FiniteVolumeField<T>::NORMAL_GRADIENT:
            case FiniteVolumeField<T>::SYMMETRY:
                break;

            default:
                throw Exception("fv", "laplacian<T>", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(cell, cell, centralCoeff);
    }

    return eqn;
}

}

#endif
