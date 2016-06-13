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

VectorFiniteVolumeField gravity(const ScalarFiniteVolumeField& rho, const Vector2D& g)
{
    VectorFiniteVolumeField gravity(rho.grid, "g");

    for(const Cell& cell: rho.grid.fluidCells())
    {
        Vector2D &gs = gravity[cell.id()] = Vector2D(0., 0.);

        for(const InteriorLink &nb: cell.neighbours())
            gs += dot(g, nb.rFaceVec())*rho.faces()[nb.face().id()]*nb.outwardNorm();

        for(const BoundaryLink &bd: cell.boundaries())
            gs += dot(g, bd.rFaceVec())*rho.faces()[bd.face().id()]*bd.outwardNorm();

        gs /= cell.volume();
    }

    interpolateFaces(gravity);

    return gravity;
}

}
