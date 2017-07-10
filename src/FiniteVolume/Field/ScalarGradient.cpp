#include "ScalarGradient.h"

ScalarGradient::ScalarGradient(const ScalarFiniteVolumeField& phi)
        :
        VectorFiniteVolumeField(phi.gridPtr(),
                                "grad" + phi.name(),
                                Vector2D(0., 0.),
                                true, false),
        phi_(phi)
{

}

void ScalarGradient::computeFaces()
{
    VectorFiniteVolumeField &gradPhi = *this;

    for(const Face& face: grid_->interiorFaces())
    {
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        gradPhi(face) = (phi_(face.rCell()) - phi_(face.lCell())) * rc / dot(rc, rc);
    }

    for(const Face& face: grid_->boundaryFaces())
    {
        Vector2D rf = face.centroid() - face.lCell().centroid();
        gradPhi(face) = (phi_(face) - phi_(face.lCell())) * rf / dot(rf, rf);
    }
}

void ScalarGradient::compute(const CellGroup& group)
{
    computeFaces();
    VectorFiniteVolumeField &gradPhi = *this;

    for (const Cell &cell: group)
    {
        Scalar sumSfx = 0., sumSfy = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D &sf = nb.outwardNorm();
            Vector2D gradF = gradPhi(nb.face());

            gradPhi(cell) += Vector2D(gradF.x * fabs(sf.x), gradF.y * fabs(sf.y));
            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D &sf = bd.outwardNorm();
            Vector2D gradF = gradPhi(bd.face());

            gradPhi(cell) += Vector2D(gradF.x * fabs(sf.x), gradF.y * fabs(sf.y));
            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        gradPhi(cell) = Vector2D(gradPhi(cell).x / sumSfx, gradPhi(cell).y / sumSfy);
    }
}

void ScalarGradient::computeWeighted(const ScalarFiniteVolumeField &w, const CellGroup &group)
{
    computeFaces();
    VectorFiniteVolumeField &gradPhi = *this;

    for (const Cell &cell: group)
    {
        Scalar sumSfx = 0., sumSfy = 0.;

        for (const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D &sf = nb.outwardNorm();
            Vector2D gradF = gradPhi(nb.face());

            gradPhi(cell) += Vector2D(gradF.x * fabs(sf.x), gradF.y * fabs(sf.y)) / w(nb.face());
            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D &sf = bd.outwardNorm();
            Vector2D gradF = gradPhi(bd.face());

            gradPhi(cell) += Vector2D(gradF.x * fabs(sf.x), gradF.y * fabs(sf.y)) / w(bd.face());
            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        gradPhi(cell) = Vector2D(gradPhi(cell).x / sumSfx, gradPhi(cell).y / sumSfy) * w(cell);
    }
}