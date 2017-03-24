#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "Equation.h"

namespace fv
{
template<typename T>
Equation<T> divexp(const VectorFiniteVolumeField &u, FiniteVolumeField<T> &field)
{
    Equation<T> eqn(field);

    for (const Cell &cell: field.grid.cellZone("fluid"))
    {
        Scalar centralCoeff = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            T faceFlux = dot(u(nb.face()), nb.outwardNorm())*field(nb.face());
            eqn.addBoundary(cell, -faceFlux);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());

            switch (field.boundaryType(bd.face()))
            {
                case FiniteVolumeField<T>::FIXED:
                    eqn.addBoundary(cell, -faceFlux * field(bd.face()));
                    break;

                case FiniteVolumeField<T>::NORMAL_GRADIENT:
                    centralCoeff += faceFlux;
                    break;

                case FiniteVolumeField<T>::SYMMETRY:
                    break;

                default:
                    throw Exception("fv", "divexp<T>", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(cell, cell, centralCoeff);
    }

    return eqn;
}

    template<typename T>
    Equation<T> div(const VectorFiniteVolumeField &u, FiniteVolumeField<T> &field)
    {
        Equation<T> eqn(field);

        for (const Cell &cell: field.grid.cellZone("fluid"))
        {
            Scalar centralCoeff = 0.;

            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar faceFlux = dot(u(nb.face()), nb.outwardNorm());

                Scalar coeff = std::min(faceFlux, 0.);
                centralCoeff += std::max(faceFlux, 0.);

                eqn.add(cell, nb.cell(), coeff);
            }

            for (const BoundaryLink &bd: cell.boundaries())
            {
                Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());

                switch (field.boundaryType(bd.face()))
                {
                    case FiniteVolumeField<T>::FIXED:
                        eqn.addBoundary(cell, -faceFlux * field(bd.face()));
                        break;

                    case FiniteVolumeField<T>::NORMAL_GRADIENT:
                        centralCoeff += faceFlux;
                        break;

                    case FiniteVolumeField<T>::SYMMETRY:
                        break;

                    default:
                        throw Exception("fv", "div<T>", "unrecognized or unspecified boundary type.");
                }
            }

            eqn.add(cell, cell, centralCoeff);
        }

        return eqn;
    }

    template<typename T>
    Equation<T> div(const ScalarFiniteVolumeField& rho, const VectorFiniteVolumeField &u, FiniteVolumeField<T> &field)
    {
        Equation<T> eqn(field);

        for (const Cell &cell: field.grid.cellZone("fluid"))
        {
            Scalar centralCoeff = 0.;

            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar faceFlux = dot(u(nb.face()), nb.outwardNorm());

                Scalar coeff = std::min(faceFlux, 0.)*rho(nb.cell());
                centralCoeff += std::max(faceFlux, 0.)*rho(cell);

                eqn.add(cell, nb.cell(), coeff);
            }

            for (const BoundaryLink &bd: cell.boundaries())
            {
                Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());

                switch (field.boundaryType(bd.face()))
                {
                    case FiniteVolumeField<T>::FIXED:
                        eqn.addBoundary(cell, -rho(bd.face()) * faceFlux * field(bd.face()));
                        break;

                    case FiniteVolumeField<T>::NORMAL_GRADIENT:
                        centralCoeff += rho(cell) * faceFlux;
                        break;

                    case FiniteVolumeField<T>::SYMMETRY:
                        break;

                    default:
                        throw Exception("fv", "div<T>", "unrecognized or unspecified boundary type.");
                }
            }

            eqn.add(cell, cell, centralCoeff);
        }

        return eqn;
    }
}

#endif
