#include "Laplacian.h"

namespace fv
{

Equation<ScalarFiniteVolumeField> laplacian(ScalarFiniteVolumeField& field)
{
    const Size nCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(5*nCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        const Index row = cell.globalIndex();
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Index col = nb.cell().globalIndex();
            const Scalar coeff = dot(nb.rCellVec(), nb.outwardNorm())/dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar coeff = dot(bd.rFaceVec(), bd.outwardNorm())/dot(bd.rFaceVec(), bd.rFaceVec());

            switch(field.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                eqn.boundaries()(row) -= coeff*field(bd.face());
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT: case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("fv", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<ScalarFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& gamma, ScalarFiniteVolumeField& field)
{
    const Size nCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(5*nCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        const Index row = cell.globalIndex();
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Index col = nb.cell().globalIndex();
            const Scalar coeff = gamma(nb.face())*dot(nb.rCellVec(), nb.outwardNorm())/dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar coeff = gamma(bd.face())*dot(bd.rFaceVec(), bd.outwardNorm())/dot(bd.rFaceVec(), bd.rFaceVec());

            switch(field.boundaryType(bd.face()))
            {
            case ScalarFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                eqn.boundaries()(row) -= coeff*field(bd.face());
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT: case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("fv", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& gamma, VectorFiniteVolumeField& field)
{
    const Size nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<VectorFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);

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

            const Scalar coeff = gamma(nb.face())*dot(nb.rCellVec(), nb.outwardNorm())/dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;

            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, colX, coeff));
            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, colY, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Scalar coeff = gamma(bd.face())*dot(bd.rFaceVec(), bd.outwardNorm())/dot(bd.rFaceVec(), bd.rFaceVec());

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
                throw Exception("fv", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, centralCoeff));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

}
