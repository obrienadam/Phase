#include "SourceEvaluation.h"

#include "FaceInterpolation.h"

namespace fv
{

VectorFiniteVolumeField source(VectorFiniteVolumeField field)
{
    for(const Cell &cell: field.grid.cellZone("fluid"))
        field(cell) *= cell.volume();

    return field;
}

VectorFiniteVolumeField inverseWeightedSource(const ScalarFiniteVolumeField& w, const VectorFiniteVolumeField& field)
{
    VectorFiniteVolumeField fb(field.grid, "fb");

    for(const Cell& cell: fb.grid.cellZone("fluid"))
    {
        Scalar sumSfx = 0., sumSfy = 0.;
        fb(cell) = Vector2D(0., 0.);

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D& sf = nb.outwardNorm();

            fb(cell) += Vector2D(field(nb.face()).x*fabs(sf.x), field(nb.face()).y*fabs(sf.y))/w(nb.face());

            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D& sf = bd.outwardNorm();

            fb(cell) += Vector2D(field(bd.face()).x*fabs(sf.x), field(bd.face()).y*fabs(sf.y))/w(bd.face());

            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        fb(cell) = w(cell)*Vector2D(fb(cell).x/sumSfx, fb(cell).y/sumSfy);
    }

    return fb;
}

VectorFiniteVolumeField gravity(const ScalarFiniteVolumeField& rho, const Vector2D& g)
{
    VectorFiniteVolumeField gravity(rho.grid, "g");
    gravity.fill(Vector2D(0., 0.));

    for(const Cell& cell: rho.grid.cellZone("fluid"))
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

}
