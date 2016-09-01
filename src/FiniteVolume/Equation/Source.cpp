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
    gravity.fill(Vector2D(0., 0.));

    //    case POTENTIAL:
    //        for(const Cell& cell: rho.grid.fluidCells())
    //        {
    //            for(const InteriorLink &nb: cell.neighbours())
    //                gravity(cell) += dot(g, nb.rFaceVec())*rho(nb.face())*nb.outwardNorm();

    //            for(const BoundaryLink &bd: cell.boundaries())
    //                gravity(cell) += dot(g, bd.rFaceVec())*rho(bd.face())*bd.outwardNorm();

    //            gravity(cell) /= cell.volume();
    //        }

    //        interpolateFaces(fv::INVERSE_VOLUME, gravity);

    for(const Cell& cell: rho.grid.fluidCells())
    {
        Scalar sumSfx = 0., sumSfy = 0.;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D& sf = nb.outwardNorm();

            gravity(cell) += rho(nb.face())*Vector2D(g.x*fabs(sf.x), g.y*fabs(sf.y));

            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D& sf = bd.outwardNorm();

            gravity(cell) += rho(bd.face())*Vector2D(g.x*fabs(sf.x), g.y*fabs(sf.y));

            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        gravity(cell) = Vector2D(gravity(cell).x/sumSfx, gravity(cell).y/sumSfy);
    }

    for(const Face& face: gravity.grid.faces())
        gravity(face) = rho(face)*g;

    return gravity;
}

ScalarFiniteVolumeField hydroStaticPressureBoundaries(const ScalarFiniteVolumeField& p, const ScalarFiniteVolumeField& rho, const Vector2D& g)
{
    //- This function corrects the pressure along boundary faces to account for the presence of gravity

    ScalarFiniteVolumeField rgh(p.grid, "rgh");
    rgh.fill(0.);

    for(const Face &face: p.grid.boundaryFaces())
    {
        if(!p.boundaryType(face) == ScalarFiniteVolumeField::NORMAL_GRADIENT)
            continue;

        const Vector2D rf = face.centroid() - face.lCell().centroid();
        const Vector2D sf = face.outwardNorm(face.lCell().centroid());

        rgh(face.lCell()) -= rho(face.lCell())*dot(g, rf)*dot(rf, sf)/dot(rf, rf);
    }

    return rgh;
}

}
