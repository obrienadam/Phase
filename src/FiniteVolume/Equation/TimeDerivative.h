#ifndef TIME_DERIVATIVE_H
#define TIME_DERIVATIVE_H

#include "Equation.h"

namespace fv
{

    template<typename T>
    Equation<T> ddt(const ScalarFiniteVolumeField &rho, FiniteVolumeField<T> &field, Scalar timeStep)
    {
        const ScalarFiniteVolumeField &rho0 = rho.oldField(0);

        Equation<T> eqn(field);

        for (const Cell &cell: field.grid.cellZone("fluid"))
        {
            eqn.add(cell, cell, rho(cell) * cell.volume() / timeStep);
            eqn.addSource(cell, -rho0(cell) * cell.volume() * field(cell) / timeStep);
        }

        return eqn;
    }

    template<typename T>
    Equation<T> ddt(FiniteVolumeField<T> &field, Scalar timeStep)
    {
        Equation<T> eqn(field);

        for (const Cell &cell: field.grid.cellZone("fluid"))
        {
            eqn.add(cell, cell, cell.volume() / timeStep);
            eqn.addSource(cell, -cell.volume() * field(cell) / timeStep);
        }

        return eqn;
    }

}

#endif
