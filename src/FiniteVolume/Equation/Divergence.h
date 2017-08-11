#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "Equation.h"
#include "JacobianField.h"

namespace fv
{
    template<typename T>
    Equation<T> div(const VectorFiniteVolumeField &u, FiniteVolumeField<T> &phi, Scalar theta = 1.)
    {
        Equation<T> eqn(phi);

        const VectorFiniteVolumeField& u0 = u.oldField(0);

        for (const Cell &cell: phi.grid().cellZone("fluid"))
        {
            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar flux = dot(u(nb.face()), nb.outwardNorm());
                Scalar flux0 = dot(u0(nb.face()), nb.outwardNorm());

                eqn.add(cell, cell, theta * std::max(flux, 0.));
                eqn.add(cell, nb.cell(), theta * std::min(flux, 0.));
                eqn.addSource(cell, (1. - theta) * std::max(flux0, 0.) * phi(cell));
                eqn.addSource(cell, (1. - theta) * std::min(flux0, 0.) * phi(nb.cell()));
            }

            for (const BoundaryLink &bd: cell.boundaries())
            {
                Scalar flux = dot(u(bd.face()), bd.outwardNorm());
                Scalar flux0 = dot(u0(bd.face()), bd.outwardNorm());

                switch (phi.boundaryType(bd.face()))
                {
                    case FiniteVolumeField<T>::FIXED:
                        eqn.addSource(cell, theta * flux * phi(bd.face()));
                        eqn.addSource(cell, (1. - theta) * flux0 * phi(bd.face()));
                        break;

                    case FiniteVolumeField<T>::NORMAL_GRADIENT:
                        eqn.add(cell, cell, theta * flux);
                        eqn.add(cell, cell, (1. - theta) * flux);
                        break;

                    case FiniteVolumeField<T>::SYMMETRY:
                        break;

                    default:
                        throw Exception("fv", "div<T>", "unrecognized or unspecified boundary type.");
                }
            }
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
