#include "CrankNicolson.h"

namespace cn
{

Equation<VectorFiniteVolumeField> div(const VectorFiniteVolumeField& u, VectorFiniteVolumeField& field)
{
    const size_t nActiveCells = field.grid.nActiveCells();

    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
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

            Scalar faceFlux = dot(u.faces()[nb.face().id()], nb.outwardNorm());

            Scalar coeff = std::min(faceFlux, 0.);
            centralCoeff += std::max(faceFlux, 0.);

            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, colX, coeff/2.));
            entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, colY, coeff/2.));

            eqn.boundaries()(rowX) -= coeff/2.*field0[nb.cell().id()].x;
            eqn.boundaries()(rowY) -= coeff/2.*field0[nb.cell().id()].y;
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar faceFlux = dot(u.faces()[bd.face().id()], bd.outwardNorm()), boundaryCoeff;

            switch(field.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.boundaries()(rowX) -= faceFlux*field.faces()[bd.face().id()].x;
                eqn.boundaries()(rowY) -= faceFlux*field.faces()[bd.face().id()].y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT: case VectorFiniteVolumeField::OUTFLOW:
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

        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowX, rowX, centralCoeff/2.));
        entries.push_back(Equation<VectorFiniteVolumeField>::Triplet(rowY, rowY, centralCoeff/2.));
        eqn.boundaries()(rowX) -= centralCoeff/2.*field0[cell.id()].x;
        eqn.boundaries()(rowY) -= centralCoeff/2.*field0[cell.id()].y;
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

}
