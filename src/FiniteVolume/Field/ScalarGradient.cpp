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

    std::fill(gradPhi.begin(), gradPhi.end(), Vector2D(0., 0.));

    for (const Cell &cell: group)
    {
        Vector2D sumSf(0., 0.);

        for (const InteriorLink &nb: cell.neighbours())
        {
            Vector2D sf = nb.outwardNorm().abs();
            gradPhi(cell) += pointwise(gradPhi(nb.face()), sf);
            sumSf += sf;
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            Vector2D sf = bd.outwardNorm().abs();
            gradPhi(cell) += pointwise(gradPhi(bd.face()), sf);
            sumSf += sf;
        }

        gradPhi(cell) = Vector2D(gradPhi(cell).x / sumSf.x, gradPhi(cell).y / sumSf.y);
    }
}