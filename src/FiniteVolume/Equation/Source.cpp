#include "Source.h"
#include "FaceInterpolation.h"

namespace fv
{

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

    interpolateFaces(fv::INVERSE_VOLUME, gravity);

    return gravity;
}

ScalarFiniteVolumeField hydroStaticPressureBoundaries(const ScalarFiniteVolumeField& p, const ScalarFiniteVolumeField& rho, const Vector2D& g)
{
    ScalarFiniteVolumeField rgh(p.grid, "rgh");
    rgh.fill(0.);

    for(const Face &face: p.grid.boundaryFaces())
    {
        if(!p.boundaryType(face.id()) == ScalarFiniteVolumeField::NORMAL_GRADIENT)
            continue;

        Vector2D rf = face.centroid() - face.lCell().centroid();
        Vector2D sf = face.outwardNorm(face.lCell().centroid());

        rgh[face.lCell().id()] -= rho[face.lCell().id()]*dot(g, rf)*dot(rf, sf)/dot(rf, rf);
    }

    return rgh;
}

}
