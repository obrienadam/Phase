#include "AdamsBashforth.h"

namespace ab
{

Equation<VectorFiniteVolumeField> laplacian(const ScalarFiniteVolumeField& gamma, VectorFiniteVolumeField& field)
{
    const size_t nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<VectorFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);
    const VectorFiniteVolumeField &field0 = field.prev(0);

    entries.reserve(10*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;

        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            size_t colX = nb.cell().globalIndex();
            size_t colY = colX + nActiveCells;

            Scalar coeff = gamma.faces()[nb.face().id()]*dot(nb.rCellVec(), nb.outwardNorm())/dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;

            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, colX, 1.5*coeff));
            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, colY, 1.5*coeff));
            eqn.boundaries()(rowX) -= -0.5*coeff*field0[nb.cell().id()].x;
            eqn.boundaries()(rowY) -= -0.5*coeff*field0[nb.cell().id()].y;
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar coeff = gamma.faces()[bd.face().id()]*dot(bd.rFaceVec(), bd.outwardNorm())/dot(bd.rFaceVec(), bd.rFaceVec());
            Vector2D source, nWall;

            switch(field.boundaryType(bd.face().id()))
            {
            case VectorFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                eqn.boundaries()(rowX) -= coeff*field.faces()[bd.face().id()].x;
                eqn.boundaries()(rowY) -= coeff*field.faces()[bd.face().id()].y;

                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT: case VectorFiniteVolumeField::OUTFLOW:
                source = bd.outwardNorm().mag()/coeff*field.boundaryRefValue(bd.face().id());
                eqn.boundaries()(rowX) -= source.x;
                eqn.boundaries()(rowY) -= source.y;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                nWall = bd.outwardNorm().unitVec();

                entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, -coeff*nWall.x*nWall.x));
                entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(rowX, rowY, -coeff*nWall.y*nWall.x));

                entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, -coeff*nWall.y*nWall.y));
                entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(rowY, rowX, -coeff*nWall.x*nWall.y));
                break;

            default:
                throw Exception("ab", "laplacian", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, 1.5*centralCoeff));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, 1.5*centralCoeff));
        eqn.boundaries()(rowX) -= -0.5*centralCoeff*field0[cell.id()].x;
        eqn.boundaries()(rowY) -= -0.5*centralCoeff*field0[cell.id()].y;
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

}
