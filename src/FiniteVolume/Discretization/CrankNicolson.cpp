#include "CrankNicolson.h"

namespace cn
{

Equation<VectorFiniteVolumeField> div(const ScalarFiniteVolumeField& rho, const VectorFiniteVolumeField& u, VectorFiniteVolumeField& field, Scalar theta)
{
    const Size nActiveCells = field.grid.nActiveCells();

    Equation<VectorFiniteVolumeField> eqn(field);
    const VectorFiniteVolumeField &field0 = field.prev(0);
    const ScalarFiniteVolumeField &rho0 = rho.prev(0);
    const VectorFiniteVolumeField &u0 = u.prev(0);

    for(const Cell& cell: field.grid.fluidCells())
    {
        const Index rowX = cell.globalIndex();
        const Index rowY = rowX + nActiveCells;
        Scalar centralCoeff0 = 0.;
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Index colX = nb.cell().globalIndex();
            const Index colY = colX + nActiveCells;

            const Scalar faceFlux0 = rho0(cell)*dot(u0(nb.face()), nb.outwardNorm());
            const Scalar faceFlux = rho(cell)*dot(u(nb.face()), nb.outwardNorm());

            Scalar coeff0 = std::min(faceFlux0, 0.);
            Scalar coeff = std::min(faceFlux, 0.);

            centralCoeff0 += std::max(faceFlux0, 0.);
            centralCoeff += std::max(faceFlux, 0.);

            eqn.add(rowX, colX, theta*coeff);
            eqn.add(rowY, colY, theta*coeff);

            eqn.boundaries()(rowX) -= (1. - theta)*coeff0*field0(nb.cell()).x;
            eqn.boundaries()(rowY) -= (1. - theta)*coeff0*field0(nb.cell()).y;
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar faceFlux = rho(cell)*dot(u(bd.face()), bd.outwardNorm());
            const Scalar faceFlux0 = rho0(cell)*dot(u0(bd.face()), bd.outwardNorm());

            switch(field.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.boundaries()(rowX) -= theta*faceFlux*field(bd.face()).x;
                eqn.boundaries()(rowY) -= theta*faceFlux*field(bd.face()).y;
                eqn.boundaries()(rowX) -= (1. - theta)*faceFlux0*field0(bd.face()).x;
                eqn.boundaries()(rowY) -= (1. - theta)*faceFlux0*field0(bd.face()).y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += faceFlux;
                centralCoeff0 += faceFlux0;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                Scalar boundaryCoeff = dot(bd.rFaceVec(), bd.outwardNorm().unitVec())/dot(bd.rFaceVec(), bd.rFaceVec());
                boundaryCoeff = 1./(1./boundaryCoeff + 1.);
                centralCoeff += faceFlux*boundaryCoeff;
            }
                break;

            default:
                throw Exception("cn", "div", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(rowX, rowX, theta*centralCoeff);
        eqn.add(rowY, rowY, theta*centralCoeff);
        eqn.boundaries()(rowX) -= (1. - theta)*centralCoeff0*field0(cell).x;
        eqn.boundaries()(rowY) -= (1. - theta)*centralCoeff0*field0(cell).y;
    }

    return eqn;
}

Equation<VectorFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& gamma, VectorFiniteVolumeField& field, Scalar theta)
{
    const Size nActiveCells = field.grid.nActiveCells();

    Equation<VectorFiniteVolumeField> eqn(field);
    const VectorFiniteVolumeField &field0 = field.prev(0);
    const ScalarFiniteVolumeField &gamma0 = gamma.prev(0);

    for(const Cell& cell: field.grid.fluidCells())
    {
        const Index rowX = cell.globalIndex();
        const Index rowY = rowX + nActiveCells;

        Scalar centralCoeff0 = 0.;
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Index colX = nb.cell().globalIndex();
            const Index colY = colX + nActiveCells;

            const Vector2D& rc = nb.rCellVec();
            const Vector2D& sf = nb.outwardNorm();

            const Scalar coeff0 = gamma0(nb.face())*dot(rc, sf)/dot(rc, rc);
            const Scalar coeff = gamma(nb.face())*dot(rc, sf)/dot(rc, rc);

            centralCoeff0 -= coeff0;
            centralCoeff -= coeff;

            eqn.add(rowX, colX, theta*coeff);
            eqn.add(rowY, colY, theta*coeff);
            eqn.boundaries()(rowX) -= (1. - theta)*coeff0*field0(nb.cell()).x;
            eqn.boundaries()(rowY) -= (1. - theta)*coeff0*field0(nb.cell()).y;
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D &rf = bd.rFaceVec();
            const Vector2D &sf = bd.outwardNorm();

            const Scalar coeff = gamma(bd.face())*dot(rf, sf)/dot(rf, rf);
            const Scalar coeff0 = gamma0(bd.face())*dot(rf, sf)/dot(rf, rf);

            switch(field.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                centralCoeff0 -= coeff;
                eqn.boundaries()(rowX) -= theta*coeff*field(bd.face()).x;
                eqn.boundaries()(rowY) -= theta*coeff*field(bd.face()).y;
                eqn.boundaries()(rowX) -= (1. - theta)*coeff0*field0(bd.face()).x;
                eqn.boundaries()(rowY) -= (1. - theta)*coeff0*field0(bd.face()).y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = bd.outwardNorm().unitVec();

                eqn.add(rowX, rowX, -coeff*nWall.x*nWall.x);
                eqn.add(rowX, rowY, -coeff*nWall.y*nWall.x);

                eqn.add(rowY, rowY, -coeff*nWall.y*nWall.y);
                eqn.add(rowY, rowX, -coeff*nWall.x*nWall.y);
            }
                break;

            default:
                throw Exception("cn", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(rowX, rowX, theta*centralCoeff);
        eqn.add(rowY, rowY, theta*centralCoeff);
        eqn.boundaries()(rowX) -= (1. - theta)*centralCoeff0*field0(cell).x;
        eqn.boundaries()(rowY) -= (1. - theta)*centralCoeff0*field0(cell).y;
    }

    return eqn;
}

}
