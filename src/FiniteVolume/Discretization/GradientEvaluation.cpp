#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

namespace fv
{

//- Scalar gradients
void
computeGradient(GradientEvaluationMethod method, const CellGroup &group, ScalarFiniteVolumeField &phi, VectorFiniteVolumeField &gradPhi)
{
    gradPhi.fill(Vector2D(0., 0.));

    for (const Face &face: phi.grid.interiorFaces())
    {
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        gradPhi(face) = (phi(face.rCell()) - phi(face.lCell())) * rc / dot(rc, rc);
    }

    for (const Face &face: phi.grid.boundaryFaces())
    {
        Vector2D rf = face.centroid() - face.lCell().centroid();
        gradPhi(face) = (phi(face) - phi(face.lCell())) * rf / dot(rf, rf);
    }

    switch (method)
    {
    case GREEN_GAUSS_CELL_CENTERED:
        for (const Cell &cell: group)
        {
            for (const InteriorLink &nb: cell.neighbours())
                gradPhi(cell) += phi(nb.face()) * nb.outwardNorm();

            for (const BoundaryLink &bd: cell.boundaries())
                gradPhi(cell) += phi(bd.face()) * bd.outwardNorm();

            gradPhi(cell) /= cell.volume();
        }

        break;

    case FACE_TO_CELL:

        for (const Cell &cell: group)
        {
            Scalar sumSfx = 0., sumSfy = 0.;
            Vector2D grad(0., 0.);

            for (const InteriorLink &nb: cell.neighbours())
            {
                const Vector2D &sf = nb.outwardNorm();

                grad = gradPhi(nb.face());

                gradPhi(cell) += Vector2D(grad.x * fabs(sf.x), grad.y * fabs(sf.y));
                sumSfx += fabs(sf.x);
                sumSfy += fabs(sf.y);
            }

            for (const BoundaryLink &bd: cell.boundaries())
            {
                const Vector2D &sf = bd.outwardNorm();

                grad = gradPhi(bd.face());

                gradPhi(cell) += Vector2D(grad.x * fabs(sf.x), grad.y * fabs(sf.y));
                sumSfx += fabs(sf.x);
                sumSfy += fabs(sf.y);
            }

            gradPhi(cell) = Vector2D(gradPhi(cell).x / sumSfx, gradPhi(cell).y / sumSfy);
        }

        break;

    default:
        throw Exception("fv", "computeGradient", "unrecognized gradient evaluation method.");
    }
}

void computeInverseWeightedGradient(const ScalarFiniteVolumeField &w, ScalarFiniteVolumeField &field,
                                    VectorFiniteVolumeField &gradField)
{
    gradField.fill(Vector2D(0., 0.));
    field.setBoundaryFaces();

    for (const Cell &cell: field.grid.cells())
    {
        Scalar sumSfx = 0., sumSfy = 0.;
        Vector2D grad(0., 0.);

        for (const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D &rc = nb.rCellVec();
            const Vector2D &sf = nb.outwardNorm();

            grad = (field(nb.cell()) - field(cell)) * rc / dot(rc, rc);

            gradField(cell) += Vector2D(grad.x * fabs(sf.x), grad.y * fabs(sf.y))/w(nb.face());
            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        for (const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D &rf = bd.rFaceVec();
            const Vector2D &sf = bd.outwardNorm();

            grad = (field(bd.face()) - field(cell)) * rf / dot(rf, rf);

            gradField(cell) += Vector2D(grad.x * fabs(sf.x), grad.y * fabs(sf.y))/w(bd.face());
            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        gradField(cell) = w(cell)*Vector2D(gradField(cell).x / sumSfx, gradField(cell).y / sumSfy);
    }

    for (const Face &face: field.grid.interiorFaces())
    {
        const Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        gradField(face) = (field(face.rCell()) - field(face.lCell())) * rc / dot(rc, rc);
    }

    for (const Face &face: field.grid.boundaryFaces())
    {
        const Vector2D rf = face.centroid() - face.lCell().centroid();
        gradField(face) = (field(face) - field(face.lCell())) * rf / dot(rf, rf);
    }
}

}
