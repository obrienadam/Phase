#include "CrankNicolson.h"

namespace cn
{

Equation<VectorFiniteVolumeField> div(const VectorFiniteVolumeField& u, VectorFiniteVolumeField& field, Scalar theta)
{
    const Size nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);
    const VectorFiniteVolumeField &field0 = field.prev(0);

    entries.reserve(10*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        const Index rowX = cell.globalIndex();
        const Index rowY = rowX + nActiveCells;
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Index colX = nb.cell().globalIndex();
            const Index colY = colX + nActiveCells;

            const Scalar faceFlux = dot(u(nb.face()), nb.outwardNorm());

            Scalar coeff = std::min(faceFlux, 0.);
            centralCoeff += std::max(faceFlux, 0.);

            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, colX, theta*coeff));
            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, colY, theta*coeff));

            eqn.boundaries()(rowX) -= (1. - theta)*coeff*field0(nb.cell()).x;
            eqn.boundaries()(rowY) -= (1. - theta)*coeff*field0(nb.cell()).y;
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());

            switch(field.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.boundaries()(rowX) -= faceFlux*field(bd.face()).x;
                eqn.boundaries()(rowY) -= faceFlux*field(bd.face()).y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += faceFlux;
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

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, theta*centralCoeff));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, theta*centralCoeff));
        eqn.boundaries()(rowX) -= (1. - theta)*centralCoeff*field0(cell).x;
        eqn.boundaries()(rowY) -= (1. - theta)*centralCoeff*field0(cell).y;
    }

    eqn.assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& gamma, VectorFiniteVolumeField& field, Scalar theta)
{
    const Size nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<VectorFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);
    const VectorFiniteVolumeField &field0 = field.prev(0);

    entries.reserve(10*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        const Index rowX = cell.globalIndex();
        const Index rowY = rowX + nActiveCells;

        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Index colX = nb.cell().globalIndex();
            const Index colY = colX + nActiveCells;

            const Vector2D& rc = nb.rCellVec();
            const Vector2D& sf = nb.outwardNorm();

            const Scalar coeff = gamma(nb.face())*dot(rc, sf)/dot(rc, rc);
            centralCoeff -= coeff;

            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, colX, theta*coeff));
            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, colY, theta*coeff));
            eqn.boundaries()(rowX) -= (1. - theta)*coeff*field0(nb.cell()).x;
            eqn.boundaries()(rowY) -= (1. - theta)*coeff*field0(nb.cell()).y;
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D &rf = bd.rFaceVec();
            const Vector2D &sf = bd.outwardNorm();

            const Scalar coeff = gamma(bd.face())*dot(rf, sf)/dot(rf, rf);

            switch(field.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                eqn.boundaries()(rowX) -= coeff*field(bd.face()).x;
                eqn.boundaries()(rowY) -= coeff*field(bd.face()).y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = bd.outwardNorm().unitVec();

                entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, -coeff*nWall.x*nWall.x));
                entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(rowX, rowY, -coeff*nWall.y*nWall.x));

                entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, -coeff*nWall.y*nWall.y));
                entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(rowY, rowX, -coeff*nWall.x*nWall.y));
            }
                break;

            default:
                throw Exception("cn", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, theta*centralCoeff));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, theta*centralCoeff));
        eqn.boundaries()(rowX) -= (1. - theta)*centralCoeff*field0(cell).x;
        eqn.boundaries()(rowY) -= (1. - theta)*centralCoeff*field0(cell).y;
    }

    eqn.assemble(entries);
    return eqn;
}

}
