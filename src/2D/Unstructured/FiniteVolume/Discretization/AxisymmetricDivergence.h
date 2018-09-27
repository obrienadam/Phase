#ifndef PHASE_AXISYMMETRIC_DIVERGENCE_H
#define PHASE_AXISYMMETRIC_DIVERGENCE_H

#include "FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace axi
{
template<class T>
FiniteVolumeEquation<T> div(const VectorFiniteVolumeField &u,
                            FiniteVolumeField<T> &phi,
                            Scalar theta = 1.)
{
    FiniteVolumeEquation<T> eqn(phi);
    const VectorFiniteVolumeField &u0 = u.oldField(0);
    const FiniteVolumeField<T> &phi0 = phi.oldField(0);

    for (const Cell &cell: phi.cells())
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.polarOutwardNorm();

            Scalar flux = theta * dot(u(nb.face()), sf);
            Scalar flux0 = (1. - theta) * dot(u0(nb.face()), sf);

            eqn.add(cell, cell, std::max(flux, 0.));
            eqn.add(cell, nb.cell(), std::min(flux, 0.));
            eqn.addSource(cell, std::max(flux0, 0.) * phi0(cell));
            eqn.addSource(cell, std::min(flux0, 0.) * phi0(nb.cell()));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.polarOutwardNorm();
            Scalar flux = theta * dot(u(bd.face()), sf);
            Scalar flux0 = (1. - theta) * dot(u0(bd.face()), sf);

            switch (phi.boundaryType(bd.face()))
            {
            case FiniteVolumeField<T>::FIXED:
                eqn.addSource(cell, flux * phi(bd.face()));
                eqn.addSource(cell, flux0 * phi0(bd.face()));
                break;

            case FiniteVolumeField<T>::NORMAL_GRADIENT:
                eqn.add(cell, cell, flux);
                eqn.addSource(cell, flux0 * phi0(bd.face()));
                break;

            case FiniteVolumeField<T>::SYMMETRY:
                break;

            default:
                throw Exception("axi", "div<T>", "unrecognized or unspecified boundary type.");
            }
        }
    }

    return eqn;
}
}

#endif
