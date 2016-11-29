#include "Divergence.h"

namespace fv
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField& u, ScalarFiniteVolumeField& field)
{
    Equation<ScalarFiniteVolumeField> eqn(field);

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

            eqn.add(row, col, coeff);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());

            switch(field.boundaryType(bd.face()))
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

        eqn.add(row, row, centralCoeff);
    }

    return eqn;
}

Equation<VectorFiniteVolumeField> div(const VectorFiniteVolumeField& u, VectorFiniteVolumeField& field)
{
    const Size nActiveCells = field.grid.nActiveCells();

    Equation<VectorFiniteVolumeField> eqn(field);

    for(const Cell& cell: field.grid.fluidCells())
    {
        Index rowX = cell.globalIndex();
        Index rowY = rowX + nActiveCells;
        Scalar centralCoeffX = 0.;
        Scalar centralCoeffY = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            Index colX = nb.cell().globalIndex();
            Index colY = colX + nActiveCells;

            Scalar faceFlux = dot(u.faces()[nb.face().id()], nb.outwardNorm());

            Scalar coeffX = std::min(faceFlux, 0.);
            Scalar coeffY = coeffX;

            centralCoeffX += std::max(faceFlux, 0.);
            centralCoeffY += std::max(faceFlux, 0.);

            eqn.add(rowX, colX, coeffX);
            eqn.add(rowY, colY, coeffY);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar faceFlux = dot(u(bd.face()), bd.outwardNorm());

            switch(field.boundaryType(bd.face()))
            {
            case VectorFiniteVolumeField::FIXED:
                eqn.boundaries()(rowX) -= faceFlux*field(bd.face()).x;
                eqn.boundaries()(rowY) -= faceFlux*field(bd.face()).y;
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeffX += faceFlux;
                centralCoeffY += faceFlux;
                break;

            case VectorFiniteVolumeField::SYMMETRY:
                break;

            default:
                throw Exception("fv", "div", "unrecognized or unspecified boundary type.");
            }
        }

        eqn.add(rowX, rowX, centralCoeffX);
        eqn.add(rowY, rowY, centralCoeffY);
    }

    return eqn;
}

}
