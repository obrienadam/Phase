#include <algorithm>

#include "AdvectionTerm.h"
#include "Exception.h"

AdvectionTerm::AdvectionTerm(const VectorFiniteVolumeField &u, const ScalarFiniteVolumeField &var)
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
            const Vector2D& uf = u.faces()[nb.face().id()];

            Scalar coeff = std::min(dot(uf, nb.outwardNorm()), 0.);
            centralCoeff += std::max(dot(uf, nb.outwardNorm()), 0.);

            coefficients_.push_back(Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            switch(var.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                sources_[row] -= dot(u.faces()[bd.face().id()], bd.outwardNorm());
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                break;

            default:
                throw Exception("AdvectionTerm", "AdvectionTerm", "unrecognized or unspecified boundary type.");
            }
        }

        coefficients_.push_back(Triplet(row, row, centralCoeff));
    }
}

AdvectionTerm::AdvectionTerm(const VectorFiniteVolumeField &u, const VectorFiniteVolumeField &var)
    :
      Term(var)
{

}

//- External functions

AdvectionTerm div(const VectorFiniteVolumeField &u, const ScalarFiniteVolumeField &var)
{
    return AdvectionTerm(u, var);
}

AdvectionTerm div(const VectorFiniteVolumeField& u, const VectorFiniteVolumeField& var)
{
    return AdvectionTerm(u, var);
}
