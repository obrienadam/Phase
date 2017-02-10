#ifndef GHOST_CELL_IMMERSED_BOUNDARY_H
#define GHOST_CELL_IMMERSED_BOUNDARY_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"

namespace gc
{

template<class ImmersedBoundaryObject_t, typename T> // Allow for use of reference wrappers as well
Equation<T> ib(const std::vector<ImmersedBoundaryObject_t>& ibObjs, FiniteVolumeField<T>& field)
{
    Equation<T> eqn(field);

    for(const ImmersedBoundaryObject& ibObj: ibObjs)
        for(const Cell &cell: ibObj.ibCells())
        {
            const Point2D& imagePoint = ibObj.imagePoint(cell);
            const std::vector< Ref<const Cell> > &kNN = ibObj.imagePointCells(cell);
            const Interpolation &bi = ibObj.imagePointInterpolation(cell);
            std::vector<Scalar> coeffs = bi(imagePoint);

            Scalar centralCoeff, l;
            // Then here we shall do the boundary assembly!
            switch(ibObj.boundaryType(field.name()))
            {
            case ImmersedBoundaryObject::FIXED:
                centralCoeff = 1.;
                break;

            case ImmersedBoundaryObject::NORMAL_GRADIENT:
                centralCoeff = -1.;
                break;

            case ImmersedBoundaryObject::CONTACT_ANGLE:
                centralCoeff = -1;
                break;

            case ImmersedBoundaryObject::PARTIAL_SLIP:
                l = (cell.centroid() - imagePoint).mag();
                centralCoeff = 1. + 2.*ibObj.boundaryRefValue(field.name())/l;

                for(Scalar &coeff: coeffs)
                    coeff *= 1. - 2.*ibObj.boundaryRefValue(field.name())/l;

                break;

            default:
                throw Exception("gc", "ib<T>", "invalid boundary type.");
            }

            eqn.add(cell, cell, centralCoeff);

            int i = 0;
            for(const Cell& bCell: kNN)
                eqn.add(cell, bCell, coeffs[i++]);
        }

    return eqn;
}

}

#endif
