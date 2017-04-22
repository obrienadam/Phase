#ifndef GHOST_CELL_IMMERSED_BOUNDARY_H
#define GHOST_CELL_IMMERSED_BOUNDARY_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"

namespace ib
{

    template<class ImmersedBoundaryObject_t, typename T>
    // Allow for use of reference wrappers as well
    Equation<T> gc(const std::vector<ImmersedBoundaryObject_t> &ibObjs, FiniteVolumeField<T> &field)
    {
        Equation <T> eqn(field);

        for (const ImmersedBoundaryObject &ibObj: ibObjs)
            for (const GhostCellStencil &stencil: ibObj.stencils())
            {
                Scalar centralCoeff, l, lambda;
                std::vector<Scalar> coeffs = stencil.ipCoeffs();

                //- Boundary assembly
                switch (ibObj.boundaryType(field.name()))
                {
                    case ImmersedBoundaryObject::FIXED:
                        centralCoeff = 0.5;
                        for (Scalar &coeff: coeffs)
                            coeff *= 0.5;

                        eqn.addBoundary(stencil.cell(), ibObj.getBoundaryRefValue<T>(field.name()));

                        break;

                    case ImmersedBoundaryObject::NORMAL_GRADIENT:
                        centralCoeff = -1.;
                        break;

                    case ImmersedBoundaryObject::CONTACT_ANGLE:
                        centralCoeff = -1;
                        break;

                    case ImmersedBoundaryObject::PARTIAL_SLIP:
                        l = stencil.length();
                        lambda = ibObj.getBoundaryRefValue<Scalar>(field.name());
                        centralCoeff = 1. + 2. * lambda / l;

                        for (Scalar &coeff: coeffs)
                            coeff *= 1. - 2. * lambda / l;

                        break;

                    default:
                        throw Exception("ib", "gc<ImmersedBoundaryObject_t, T>", "invalid boundary type.");
                }

                eqn.add(stencil.cell(), stencil.cell(), centralCoeff);

                int i = 0;
                for (const Cell &ipCell: stencil.ipCells())
                    eqn.add(stencil.cell(), ipCell, coeffs[i++]);
            }

        return eqn;
    }

    template<class ImmersedBoundaryObject_t>
    // Allow for use of reference wrappers as well
    Equation<Scalar> cl(const std::vector<ImmersedBoundaryObject_t> &ibObjs, const VectorFiniteVolumeField& gradGamma,
                           ScalarFiniteVolumeField &gamma, Scalar theta)
    {
        Equation<Scalar> eqn(gamma);

        for (const ImmersedBoundaryObject &ibObj: ibObjs)
            for (const GhostCellStencil &stencil: ibObj.stencils())
            {
                std::vector<Scalar> coeffs = stencil.ipCoeffs();

                eqn.add(stencil.cell(), stencil.cell(), -1./stencil.length());
                eqn.addBoundary(stencil.cell(), -stencil.ipValue(gradGamma).mag()*cos(theta));

                int i = 0;
                for (const Cell &ipCell: stencil.ipCells())
                    eqn.add(stencil.cell(), ipCell, coeffs[i++]/stencil.length());
            }

        return eqn;
    }

}

#endif
