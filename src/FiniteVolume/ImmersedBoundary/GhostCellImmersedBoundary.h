#ifndef GHOST_CELL_IMMERSED_BOUNDARY_H
#define GHOST_CELL_IMMERSED_BOUNDARY_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"
#include "ForcingCellStencil.h"

namespace ib
{

    template<class ImmersedBoundaryObject_t, typename T>
    // Allow for use of reference wrappers as well
    Equation<T> gc(const std::vector<ImmersedBoundaryObject_t> &ibObjs, FiniteVolumeField<T> &field)
    {
        Equation <T> eqn(field);

//        for (const ImmersedBoundaryObject &ibObj: ibObjs)

//            for (const GhostCellStencil &stencil: ibObj.ghostCellStencils())
//            {
//                Scalar centralCoeff, l, lambda;
//                std::vector<Scalar> coeffs = stencil.ipCoeffs();

//                //- Boundary assembly
//                switch (ibObj.boundaryType(field.name()))
//                {
//                    case ImmersedBoundaryObject::FIXED:
//                        centralCoeff = 0.5;
//                        for (Scalar &coeff: coeffs)
//                            coeff *= 0.5;

//                        eqn.addSource(stencil.cell(), -ibObj.getBoundaryRefValue<T>(field.name()));

//                        break;

//                    case ImmersedBoundaryObject::NORMAL_GRADIENT:
//                        centralCoeff = -1.;
//                        break;

//                    case ImmersedBoundaryObject::CONTACT_ANGLE:
//                        centralCoeff = -1;
//                        break;

//                    case ImmersedBoundaryObject::PARTIAL_SLIP:
//                        l = stencil.length();
//                        lambda = ibObj.getBoundaryRefValue<Scalar>(field.name());
//                        centralCoeff = 1. + 2. * lambda / l;

//                        for (Scalar &coeff: coeffs)
//                            coeff *= 1. - 2. * lambda / l;

//                        break;

//                    default:
//                        throw Exception("ib", "gc<ImmersedBoundaryObject_t, T>", "invalid boundary type.");
//                }

//                eqn.add(stencil.cell(), stencil.cell(), centralCoeff);

//                int i = 0;
//                for (const Cell &ipCell: stencil.ipCells())
//                    eqn.add(stencil.cell(), ipCell, coeffs[i++]);

//                for(const Cell& cell: ibObj.freshlyClearedCells())
//                {
//                    ForcingCellStencil st(cell, ibObj.shape(), field.grid.cellZone("fluid"));

//                    eqn.add(cell, cell, 1.);
//                    eqn.add(cell, st.iCells()[0], -st.iCoeffs()[0]);
//                    eqn.add(cell, st.iCells()[1], -st.iCoeffs()[1]);
//                    eqn.addSource(cell, -st.bCoeff()*ibObj.getBoundaryRefValue<T>(field.name()));
//                }
//            }

        return eqn;
    }
}

#endif
