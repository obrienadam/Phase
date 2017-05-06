#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include "Equation.h"

namespace cn
{
template<typename T>
Equation<T> div(const VectorFiniteVolumeField &u, FiniteVolumeField<T> &field,
                Scalar theta = 0.5)
{
    Equation<T> eqn(field);
    const VectorFiniteVolumeField &u0 = u.prev(0);

    for (const Cell &cell: field.grid.cellZone("fluid"))
    {
        Scalar centralCoeff0 = 0.;
        Scalar centralCoeff = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar faceFlux0 = dot(u0(nb.face()), nb.outwardNorm());
            Scalar faceFlux = dot(u(nb.face()), nb.outwardNorm());

            Scalar coeff0 = std::min(faceFlux0, 0.);
            Scalar coeff = std::min(faceFlux, 0.);

            centralCoeff0 += std::max(faceFlux0, 0.);
            centralCoeff += std::max(faceFlux, 0.);

            eqn.add(cell, nb.cell(), theta * coeff);
            eqn.addBoundary(cell, -(1. - theta) * coeff0 * field(nb.cell()));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());
            const Scalar faceFlux0 = dot(u0(bd.face()), bd.outwardNorm());

            switch (field.boundaryType(bd.face()))
            {
            case FiniteVolumeField<T>::FIXED:
                eqn.addBoundary(cell, -theta * faceFlux * field(bd.face()));
                eqn.addBoundary(cell, -(1. - theta) * faceFlux0 * field(bd.face()));
                break;

            case FiniteVolumeField<T>::NORMAL_GRADIENT:
                centralCoeff += faceFlux;
                centralCoeff0 += faceFlux0;
                break;

            case FiniteVolumeField<T>::SYMMETRY:
                break;

            default:
                throw Exception("cn", "div<T>", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(cell, cell, theta * centralCoeff);
        eqn.addBoundary(cell, -(1. - theta) * centralCoeff0 * field(cell));
    }

    return eqn;
}

template<typename T>
Equation<T> div(const ScalarFiniteVolumeField &rho, const VectorFiniteVolumeField &u, FiniteVolumeField<T> &field,
                Scalar theta = 0.5)
{
    Equation<T> eqn(field);
    const ScalarFiniteVolumeField &rho0 = rho.prev(0);
    const VectorFiniteVolumeField &u0 = u.prev(0);

    for (const Cell &cell: field.grid.cellZone("fluid"))
    {
        Scalar centralCoeff0 = 0.;
        Scalar centralCoeff = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            Scalar faceFlux0 = dot(u0(nb.face()), nb.outwardNorm());
            Scalar faceFlux = dot(u(nb.face()), nb.outwardNorm());

            Scalar coeff0 = std::min(faceFlux0, 0.)*rho0(nb.cell());
            Scalar coeff = std::min(faceFlux, 0.)*rho(nb.cell());

            centralCoeff0 += std::max(faceFlux0, 0.)*rho0(cell);
            centralCoeff += std::max(faceFlux, 0.)*rho(cell);

            eqn.add(cell, nb.cell(), theta * coeff);
            eqn.addSource(cell, (1. - theta) * coeff0 * field(nb.cell()));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());
            const Scalar faceFlux0 = dot(u0(bd.face()), bd.outwardNorm());

            switch (field.boundaryType(bd.face()))
            {
            case FiniteVolumeField<T>::FIXED:
                eqn.addSource(cell, theta * faceFlux * rho(bd.face()) * field(bd.face()));
                eqn.addSource(cell, (1. - theta) * faceFlux0 * rho0(bd.face()) * field(bd.face()));
                break;

            case FiniteVolumeField<T>::NORMAL_GRADIENT:
                centralCoeff += faceFlux * rho(cell);
                centralCoeff0 += faceFlux0 * rho0(cell);
                break;

            case FiniteVolumeField<T>::SYMMETRY: //- No flux through the boundary face
                break;

            default:
                throw Exception("cn", "div<T>", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(cell, cell, theta * centralCoeff);
        eqn.addSource(cell, (1. - theta) * centralCoeff0 * field(cell));
    }

    return eqn;
}

template<typename T>
Equation<T> laplacian(const ScalarFiniteVolumeField &gamma, FiniteVolumeField<T> &field, Scalar theta = 0.5)
{
    Equation<T> eqn(field);
    const ScalarFiniteVolumeField &gamma0 = gamma.prev(0);

    for (const Cell &cell: field.grid.cellZone("fluid"))
    {
        Scalar centralCoeff0 = 0.;
        Scalar centralCoeff = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D &rc = nb.rCellVec();
            const Vector2D &sf = nb.outwardNorm();

            Scalar coeff0 = gamma0(nb.face()) * dot(rc, sf) / dot(rc, rc);
            Scalar coeff = gamma(nb.face()) * dot(rc, sf) / dot(rc, rc);

            centralCoeff0 -= coeff0;
            centralCoeff -= coeff;

            eqn.add(cell, nb.cell(), theta * coeff);
            eqn.addSource(cell, (1. - theta) * coeff0 * field(nb.cell()));
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D &rf = bd.rFaceVec();
            const Vector2D &sf = bd.outwardNorm();

            const Scalar coeff = gamma(bd.face()) * dot(rf, sf) / dot(rf, rf);
            const Scalar coeff0 = gamma0(bd.face()) * dot(rf, sf) / dot(rf, rf);

            switch (field.boundaryType(bd.face()))
            {
            case FiniteVolumeField<T>::FIXED:
                centralCoeff -= coeff;
                centralCoeff0 -= coeff;
                eqn.addSource(cell, theta * coeff * field(bd.face()));
                eqn.addSource(cell, (1. - theta) * coeff0 * field(bd.face()));
                break;

            case FiniteVolumeField<T>::NORMAL_GRADIENT:
            case FiniteVolumeField<T>::SYMMETRY: //- May need to be modified for vector fields, can have some normal component of shear
                break;
            default:
                throw Exception("cn", "laplacian<T>", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(cell, cell, theta * centralCoeff);
        eqn.addSource(cell, (1. - theta) * centralCoeff0 * field(cell));
    }

    return eqn;
}

}

#endif
