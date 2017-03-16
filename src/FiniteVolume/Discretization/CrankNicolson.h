#ifndef CRANK_NICOLSON_H
#define CRANK_NICOLSON_H

#include "Equation.h"

namespace cn
{

    template<typename T>
    Equation<T> div(const ScalarFiniteVolumeField &rho, const VectorFiniteVolumeField &u, FiniteVolumeField<T> &field,
                    Scalar theta = 0.5)
    {
        Equation<T> eqn(field);
        const ScalarFiniteVolumeField &rho0 = rho.prev(0);
        const VectorFiniteVolumeField &u0 = u.prev(0);
        const FiniteVolumeField<T> &field0 = field.prev(0);

        for (const Cell &cell: field.grid.cellZone("fluid"))
        {
            Scalar centralCoeff0 = 0.;
            Scalar centralCoeff = 0.;

            for (const InteriorLink &nb: cell.neighbours())
            {
                Scalar faceFlux0 = rho0(nb.face()) * dot(u0(nb.face()), nb.outwardNorm());
                Scalar faceFlux = rho(nb.face()) * dot(u(nb.face()), nb.outwardNorm());

                Scalar coeff0 = std::min(faceFlux0, 0.);
                Scalar coeff = std::min(faceFlux, 0.);

                centralCoeff0 += std::max(faceFlux0, 0.);
                centralCoeff += std::max(faceFlux, 0.);

                eqn.add(cell, nb.cell(), theta * coeff);
                eqn.addBoundary(cell, -(1. - theta) * coeff0 * field0(nb.cell()));
            }

            for (const BoundaryLink &bd: cell.boundaries())
            {
                const Scalar faceFlux = rho(bd.face()) * dot(u(bd.face()), bd.outwardNorm());
                const Scalar faceFlux0 = rho0(bd.face()) * dot(u0(bd.face()), bd.outwardNorm());

                switch (field.boundaryType(bd.face()))
                {
                    case FiniteVolumeField<T>::FIXED:
                        eqn.addBoundary(cell, -theta * faceFlux * field(bd.face()));
                        eqn.addBoundary(cell, -(1. - theta) * faceFlux0 * field0(bd.face()));
                        break;

                    case FiniteVolumeField<T>::NORMAL_GRADIENT:
                        centralCoeff += faceFlux;
                        centralCoeff0 += faceFlux0;
                        break;

                    case FiniteVolumeField<T>::SYMMETRY:
                        eqn.addBoundary(cell, -theta * faceFlux * field(bd.face()));
                        eqn.addBoundary(cell, -(1. - theta) * faceFlux0 * field0(bd.face()));
                        break;

                    default:
                        throw Exception("cn", "div<T>", "unrecognized or unspecified boundary type.");
                }
            }

            eqn.add(cell, cell, theta * centralCoeff);
            eqn.addBoundary(cell, -(1. - theta) * centralCoeff0 * field0(cell));
        }

        return eqn;
    }

    template<typename T>
    Equation<T> laplacian(const ScalarFiniteVolumeField &gamma, FiniteVolumeField<T> &field, Scalar theta = 0.5)
    {
        Equation<T> eqn(field);
        const ScalarFiniteVolumeField &gamma0 = gamma.prev(0);
        const FiniteVolumeField<T> &field0 = field.prev(0);

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
                eqn.addBoundary(cell, -(1. - theta) * coeff0 * field0(nb.cell()));
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
                        eqn.addBoundary(cell, -theta * coeff * field(bd.face()));
                        eqn.addBoundary(cell, -(1. - theta) * coeff0 * field0(bd.face()));
                        break;

                    case FiniteVolumeField<T>::NORMAL_GRADIENT:
                    case FiniteVolumeField<T>::SYMMETRY:
                        break;
                    default:
                        throw Exception("cn", "laplacian<T>", "unrecognized or unspecified boundary type.");
                }
            }

            eqn.add(cell, cell, theta * centralCoeff);
            eqn.addBoundary(cell, -(1. - theta) * centralCoeff0 * field0(cell));
        }

        return eqn;
    }

}

#endif
