#ifndef TIME_DERIVATIVE_H
#define TIME_DERIVATIVE_H

#include "Equation.h"

namespace fv
{
    template<typename T>
    Equation<T> ddt(Scalar rho, FiniteVolumeField<T>& field, Scalar timeStep, const CellGroup& cells)
    {
        Equation<T> eqn(field);

        for (const Cell &cell: cells)
        {
            eqn.add(cell, cell, rho * cell.volume() / timeStep);
            eqn.addSource(cell, -rho * cell.volume() * field(cell) / timeStep);
        }

        return eqn;
    }

    template<typename T>
    Equation<T> ddt(const ScalarFiniteVolumeField &rho, FiniteVolumeField<T> &field, Scalar timeStep, const CellGroup& cells)
    {
        const ScalarFiniteVolumeField &rho0 = rho.oldField(0);

        Equation<T> eqn(field);

        for (const Cell &cell: cells)
        {
            eqn.add(cell, cell, rho(cell) * cell.volume() / timeStep);
            eqn.addSource(cell, -rho0(cell) * cell.volume() * field(cell) / timeStep);
        }

        return eqn;
    }

    template<typename T>
    Equation<T> ddt(FiniteVolumeField<T> &field, Scalar timeStep, const CellGroup& cells)
    {
        Equation<T> eqn(field);

        for (const Cell &cell: cells)
        {
            eqn.add(cell, cell, cell.volume() / timeStep);
            eqn.addSource(cell, -cell.volume() * field(cell) / timeStep);
        }

        return eqn;
    }

    template <class T>
    Equation<T> ddt(Scalar rho, FiniteVolumeField<T>& field, Scalar timeStep)
    {
        return ddt(rho, field, timeStep, field.grid()->cellZone("fluid"));
    }

    template <class T>
    Equation<T> ddt(const ScalarFiniteVolumeField& rho, FiniteVolumeField<T>& field, Scalar timeStep)
    {
        return ddt(rho, field, timeStep, field.grid()->cellZone("fluid"));
    }

    template <class T>
    Equation<T> ddt(FiniteVolumeField<T>& field, Scalar timeStep)
    {
        return ddt(field, timeStep, field.grid()->cellZone("fluid"));
    }
}

#endif
