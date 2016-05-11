#include "Plic.h"

namespace plic
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);
    VectorFiniteVolumeField gradField = grad(field); //- should this be smoothed?

    entries.reserve(5*field.grid.nActiveCells());

    const Scalar toler = 1e-1;
    // Figure out the plic reconstruction
    for(const Cell &cell: field.grid.fluidCells())
    {
        bool isInterfaceCell = field[cell.id()] > toler && field[cell.id()] < 1. - 1e-10;

        Scalar centralCoeff = 0.;
        const size_t row = cell.globalIndex();

        for(const InteriorLink &nb: cell.neighbours())
        {
            size_t col = nb.cell().globalIndex();
            Scalar coeff;

            if(isInterfaceCell)
            {
                Vector2D norm = -gradField[cell.id()];
                Line2D line(cell.centroid(), norm.normalVec());
            }
            else
            {
                Scalar faceFlux = dot(u.faces()[nb.face().id()], nb.outwardNorm());

                coeff = std::min(faceFlux, 0.);
                centralCoeff += std::max(faceFlux, 0.);
            }

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {

        }
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

}
