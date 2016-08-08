#include "Divergence.h"

namespace fv
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField& u, ScalarFiniteVolumeField& field)
{
    const Size nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(5*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        Index row = cell.globalIndex();
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            Index col = nb.cell().globalIndex();

            Scalar faceFlux = dot(u.faces()[nb.face().id()], nb.outwardNorm());

            Scalar coeff = std::min(faceFlux, 0.);
            centralCoeff += std::max(faceFlux, 0.);

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar faceFlux = dot(u.faces()[bd.face().id()], bd.outwardNorm());

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= faceFlux*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += faceFlux;
                break;

            case ScalarFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("fv", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(row, row, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Equation<VectorFiniteVolumeField> div(const VectorFiniteVolumeField& u, VectorFiniteVolumeField& field)
{
    const Size nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);

    entries.reserve(5*nActiveCells);

    for(const Cell& cell: field.grid.fluidCells())
    {
        Index rowX = cell.globalIndex();
        Index rowY = rowX + nActiveCells;
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            Index colX = nb.cell().globalIndex();
            Index colY = colX + nActiveCells;

            Scalar faceFlux = dot(u.faces()[nb.face().id()], nb.outwardNorm());

            Scalar coeff = std::min(faceFlux, 0.);
            centralCoeff += std::max(faceFlux, 0.);

            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, colX, coeff));
            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, colY, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar faceFlux = dot(u.faces()[bd.face().id()], bd.outwardNorm()), boundaryCoeff;

            switch(field.boundaryType(bd.face().id()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.boundaries()(rowX) -= faceFlux*field.faces()[bd.face().id()].x;
                eqn.boundaries()(rowY) -= faceFlux*field.faces()[bd.face().id()].y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += faceFlux;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                boundaryCoeff = dot(bd.rFaceVec(), bd.outwardNorm().unitVec())/dot(bd.rFaceVec(), bd.rFaceVec());
                boundaryCoeff = 1./(1./boundaryCoeff + 1.);
                centralCoeff += faceFlux*boundaryCoeff;
                break;

            default:
                throw Exception("fv", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, centralCoeff));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

}
