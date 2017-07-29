#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "Equation.h"
#include "JacobianField.h"

namespace fv
{
    template<typename T>
    Equation<T> div(const VectorFiniteVolumeField &u, FiniteVolumeField<T> &field)
    {
        Equation<T> eqn(field);

        for (const Cell &cell: field.grid().cellZone("fluid"))
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
                        eqn.addSource(cell, faceFlux * field(bd.face()));
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

    template<class T>
    Equation<T> divc(const VectorFiniteVolumeField& u, FiniteVolumeField<T>& field)
    {
        Equation<T> eqn(field);

        for (const Cell &cell: field.grid().cellZone("fluid"))
        {
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar faceFlux = dot(u(nb.face()), nb.outwardNorm());
                Scalar g = nb.cell().volume() / (cell.volume() + nb.cell().volume());
                eqn.add(cell, cell, g*faceFlux);
                eqn.add(cell, nb.cell(), (1. - g)*faceFlux);
            }

            for (const BoundaryLink &bd: cell.boundaries())
            {
                Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());

                switch (field.boundaryType(bd.face()))
                {
                    case FiniteVolumeField<T>::FIXED:
                        eqn.addSource(cell, faceFlux * field(bd.face()));
                        break;

                    case FiniteVolumeField<T>::NORMAL_GRADIENT:
                        eqn.add(cell, cell, faceFlux);
                        break;

                    case FiniteVolumeField<T>::SYMMETRY:
                        break;

                    default:
                        throw Exception("fv", "divc<T>", "unrecognized or unspecified boundary type.");
                }
            }
        }

        return eqn;
    }

    Equation<Vector2D> div(const VectorFiniteVolumeField& phiU,
                           const JacobianField& gradU,
                           VectorFiniteVolumeField &u);
}

#endif
