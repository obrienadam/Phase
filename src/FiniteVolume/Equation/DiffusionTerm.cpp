#include "DiffusionTerm.h"
#include "Exception.h"

DiffusionTerm::DiffusionTerm(const ScalarFiniteVolumeField& var)
    :
      Term(var.grid)
{
    sources_.resize(var.grid.nActiveCells(), 0.);

    for(const Cell& cell: var.grid.cells)
    {
        size_t row = cell.globalIndex();
        Scalar centralCoeff = 0.;

        const auto &neighbours = cell.neighbours();

        for(const InteriorLink &nb: neighbours)
        {
            size_t col = nb.cell().globalIndex();
            Scalar coeff = dot(nb.rCellVec(), nb.outwardNorm())/dot(nb.rCellVec(), nb.rCellVec());
            centralCoeff -= coeff;

            coefficients_.push_back(Triplet(row, col, coeff));
        }

        const auto &boundaries = cell.boundaries();

        for(const BoundaryLink &bd: boundaries)
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

//- External functions
DiffusionTerm laplacian(const ScalarFiniteVolumeField &var)
{
    return DiffusionTerm(var);
}
