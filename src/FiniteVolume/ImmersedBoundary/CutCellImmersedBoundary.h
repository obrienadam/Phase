#ifndef CUT_CELL_IMMERSED_BOUNDARY_H
#define CUT_CELL_IMMERSED_BOUNDARY_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"

namespace ib
{


template<class ImmersedBoundaryObject_t, typename T>
Equation<T> div(const std::vector<ImmersedBoundaryObject_t>& ibObjs, const VectorFiniteVolumeField& u, FiniteVolumeField<T>& field)
{
    Equation<T> eqn(field);


    return eqn;
}

}

#endif
