#include "DiffusionTerm.h"
#include "Exception.h"

DiffusionTerm::DiffusionTerm(const ScalarFiniteVolumeField& var)
    :
      Term(var)
{
    for(const Cell& cell: var.grid.cells)
    {
        size_t row = cell.globalIndex();
        Scalar centralCoeff = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            size_t col = nb.cell().globalIndex();
            Scalar coeff = dot(nb.rCellVec(), nb.outwardNorm())/dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;

            coefficients_.push_back(Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar coeff = dot(bd.rFaceVec(), bd.outwardNorm())/dot(bd.rFaceVec(), bd.rFaceVec());

            switch(var.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                sources_[row] -= coeff*var.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                break;

            default:
                throw Exception("DiffusionTerm", "DiffusionTerm", "unrecognized or unspecified boundary type.");
            }
        }

        coefficients_.push_back(Triplet(row, row, centralCoeff));
    }
}

DiffusionTerm::DiffusionTerm(const VectorFiniteVolumeField &var)
    :
      Term(var)
{
    const size_t nActiveCells = var.grid.nActiveCells();

    for(const Cell& cell: var.grid.cells)
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + nActiveCells;
        Scalar centralCoeff = 0.;

        for(const InteriorLink& nb: cell.neighbours())
        {
            size_t colX = nb.cell().globalIndex();
            size_t colY = colX + nActiveCells;

            Scalar coeff = dot(nb.rCellVec(), nb.outwardNorm())/dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;

            coefficients_.push_back(Triplet(rowX, colX, coeff));
            coefficients_.push_back(Triplet(rowY, colY, coeff));
        }

        for(const BoundaryLink& bd: cell.boundaries())
        {
            Scalar coeff = dot(bd.rFaceVec(), bd.outwardNorm())/dot(bd.rFaceVec(), bd.rFaceVec());

            switch(var.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                centralCoeff -= coeff;
                sources_[rowX] -= coeff*var.faces()[bd.face().id()].x;
                sources_[rowY] -= coeff*var.faces()[bd.face().id()].y;
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                break;

            default:
                throw Exception("DiffusionTerm", "DiffusionTerm", "unrecognized or unspecified boundary type.");
            }
        }

        coefficients_.push_back(Triplet(rowX, rowX, centralCoeff));
        coefficients_.push_back(Triplet(rowY, rowY, centralCoeff));
    }
}

//- External functions
DiffusionTerm laplacian(const ScalarFiniteVolumeField &var)
{
    return DiffusionTerm(var);
}

DiffusionTerm laplacian(const VectorFiniteVolumeField &var)
{
    return DiffusionTerm(var);
}
