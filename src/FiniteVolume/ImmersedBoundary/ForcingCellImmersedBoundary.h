#ifndef FORCING_CELL_IMMERSED_BOUNDARY_H
#define FORCING_CELL_IMMERSED_BOUNDARY_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"
#include "ForcingCellStencil.h"

namespace ib
{

template<class ImmersedBoundaryObject_t, typename T>
Equation<T> dirichlet(const std::vector<ImmersedBoundaryObject_t> &ibObjs, FiniteVolumeField<T>& field)
{
    Equation<T> eqn(field);

    for(const ImmersedBoundaryObject& ibObj: ibObjs)
    {
        for(const ForcingCellStencil& st: ibObj.forcingStencils())
        {
            T bRef = ibObj.getBoundaryRefValue<T>(field.name);

            eqn.add(st.cell(), st.cell(), 1.);
            eqn.add(st.cell(), st.iCells()[0], -st.iCoeffs()[0]);
            eqn.add(st.cell(), st.iCells()[1], -st.iCoeffs()[1]);
            eqn.addSource(st.cell(), -st.bCoeff()*bRef);
        }
    }
}

template<class ImmersedBoundaryObject_t, typename T>
Equation<T> neumann(const std::vector<ImmersedBoundaryObject_t> &ibObjs, FiniteVolumeField<T>& field)
{
    Equation<T> eqn(field);

    for(const ImmersedBoundaryObject& ibObj: ibObjs)
    {
        for(const ForcingCellStencil& st: ibObj.forcingStencils())
        {
            T bRef = ibObj.getBoundaryRefValue<T>(field.name);
        }
    }
}

}

#endif
