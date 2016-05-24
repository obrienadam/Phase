#include "Source.h"

namespace fv
{

VectorFiniteVolumeField grad(const ScalarFiniteVolumeField &field)
{
    VectorFiniteVolumeField gradField(field.grid, "grad_" + field.name);

    for(const Cell& cell: gradField.grid.fluidCells())
    {
        Vector2D &gradVec = gradField[cell.id()];

        for(const InteriorLink& nb: cell.neighbours())
            gradVec += field.faces()[nb.face().id()]*nb.outwardNorm();

        for(const BoundaryLink& bd: cell.boundaries())
            gradVec += field.faces()[bd.face().id()]*bd.outwardNorm();
    }

    return gradField;
}

VectorFiniteVolumeField source(VectorFiniteVolumeField field)
{
    for(const Cell &cell: field.grid.fluidCells())
        field[cell.id()] *= cell.volume();

    return field;
}

}
