#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

namespace fv
{

//- Scalar gradients

void computeGradient(GradientEvaluationMethod method, ScalarFiniteVolumeField& field, VectorFiniteVolumeField& gradField, bool useCurrentFaceValues)
{
    gradField.fill(Vector2D(0., 0.));

    switch(method)
    {
    case GREEN_GAUSS_CELL_CENTERED:
        if(!useCurrentFaceValues)
            interpolateFaces(INVERSE_VOLUME, field);

        for(const Cell& cell: field.grid.cellZone("fluid"))
        {
            Vector2D &grad = gradField[cell.id()];

            for(const InteriorLink& nb: cell.neighbours())
                grad += field.faces()[nb.face().id()]*nb.outwardNorm();

            for(const BoundaryLink& bd: cell.boundaries())
                grad += field.faces()[bd.face().id()]*bd.outwardNorm();

            grad /= cell.volume();
        }

        break;

    case GREEN_GAUSS_NODE_CENTERED:
        if(!useCurrentFaceValues)
            interpolateNodes(field);

        for(const Cell& cell: field.grid.cellZone("fluid"))
        {
            Vector2D &grad = gradField[cell.id()];

            for(const InteriorLink& nb: cell.neighbours())
                grad += field.faces()[nb.face().id()]*nb.outwardNorm();

            for(const BoundaryLink& bd: cell.boundaries())
                grad += field.faces()[bd.face().id()]*bd.outwardNorm();

            grad /= cell.volume();
        }

        break;

    case FACE_TO_CELL:
        if(!useCurrentFaceValues)
            interpolateFaces(fv::INVERSE_VOLUME, field);

        for(const Cell &cell: field.grid.cellZone("fluid"))
        {
            Scalar sumSfx = 0., sumSfy = 0.;
            Vector2D grad(0., 0.);

            for(const InteriorLink &nb: cell.neighbours())
            {
                const Vector2D& rc = nb.rCellVec();
                const Vector2D& sf = nb.outwardNorm();

                grad = (field(nb.cell()) - field(cell))*rc/dot(rc, rc);

                gradField(cell) += Vector2D(grad.x*fabs(sf.x), grad.y*fabs(sf.y));
                sumSfx += fabs(sf.x);
                sumSfy += fabs(sf.y);
            }

            for(const BoundaryLink &bd: cell.boundaries())
            {
                const Vector2D& rf = bd.rFaceVec();
                const Vector2D& sf = bd.outwardNorm();

                grad = (field(bd.face()) - field(cell))*rf/dot(rf, rf);

                gradField(cell) += Vector2D(grad.x*fabs(sf.x), grad.y*fabs(sf.y));
                sumSfx += fabs(sf.x);
                sumSfy += fabs(sf.y);
            }

            gradField(cell) = Vector2D(gradField(cell).x/sumSfx, gradField(cell).y/sumSfy);
        }

        for(const Face& face: field.grid.interiorFaces())
        {
            const Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
            gradField(face) = (field(face.rCell()) - field(face.lCell()))*rc/dot(rc, rc);
        }

        for(const Face& face: field.grid.boundaryFaces())
        {
            const Vector2D rf = face.centroid() - face.lCell().centroid();
            gradField(face) = (field(face) - field(face.lCell()))*rf/dot(rf, rf);
        }

        break;

    default:
        throw Exception("fv", "computeGradient", "unrecognized gradient evaluation method.");
    }
}

void computeInverseWeightedGradient(const ScalarFiniteVolumeField& w, ScalarFiniteVolumeField& field, VectorFiniteVolumeField& gradField)
{
    gradField.fill(Vector2D(0., 0.));

    for(const Cell &cell: field.grid.cellZone("fluid"))
    {
        Scalar sumSfx = 0., sumSfy = 0.;
        Vector2D grad(0., 0.);

        for(const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D& rc = nb.rCellVec();
            const Vector2D& sf = nb.outwardNorm();

            grad = (field(nb.cell()) - field(cell))*rc/dot(rc, rc);

            gradField(cell) += Vector2D(grad.x*fabs(sf.x), grad.y*fabs(sf.y))/w(nb.face());
            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D& rf = bd.rFaceVec();
            const Vector2D& sf = bd.outwardNorm();

            grad = (field(bd.face()) - field(cell))*rf/dot(rf, rf);

            gradField(cell) += Vector2D(grad.x*fabs(sf.x), grad.y*fabs(sf.y))/w(bd.face());
            sumSfx += fabs(sf.x);
            sumSfy += fabs(sf.y);
        }

        gradField(cell) = w(cell)*Vector2D(gradField(cell).x/sumSfx, gradField(cell).y/sumSfy);
    }

    for(const Face& face: field.grid.interiorFaces())
    {
        const Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
        gradField(face) = (field(face.rCell()) - field(face.lCell()))*rc/dot(rc, rc);
    }

    for(const Face& face: field.grid.boundaryFaces())
    {
        const Vector2D rf = face.centroid() - face.lCell().centroid();
        gradField(face) = (field(face) - field(face.lCell()))*rf/dot(rf, rf);
    }
}

//- Vector Gradients

void computeGradient(GradientEvaluationMethod method, VectorFiniteVolumeField& field, TensorFiniteVolumeField& gradField, bool useCurrentFaceValues)
{
    gradField.fill(Tensor2D(0., 0., 0., 0.));

    switch(method)
    {
    case GREEN_GAUSS_CELL_CENTERED:
        if(!useCurrentFaceValues)
            interpolateFaces(INVERSE_VOLUME, field);

        for(const Cell& cell: field.grid.cellZone("fluid"))
        {
            Tensor2D &grad = gradField(cell);

            for(const InteriorLink& nb: cell.neighbours())
                grad += outer(field(nb.face()), nb.outwardNorm());

            for(const BoundaryLink& bd: cell.boundaries())
                grad += outer(field(bd.face()), bd.outwardNorm());

            grad /= cell.volume();
        }

        break;

    case FACE_TO_CELL:
        if(!useCurrentFaceValues)
            interpolateFaces(fv::INVERSE_VOLUME, field);

        for(const Cell &cell: field.grid.cellZone("fluid"))
        {
            Scalar sumSfx = 0., sumSfy = 0.;
            Tensor2D grad(0., 0., 0., 0.);

            for(const InteriorLink &nb: cell.neighbours())
            {
                const Vector2D& rc = nb.rCellVec();
                const Vector2D& sf = nb.outwardNorm();

                grad = outer(field(nb.cell()) - field(cell), rc/dot(rc, rc));

                gradField(cell) += Tensor2D(grad.xx*fabs(sf.x), grad.xy*fabs(sf.y),
                                            grad.yx*fabs(sf.x), grad.yy*fabs(sf.y));
                sumSfx += fabs(sf.x);
                sumSfy += fabs(sf.y);
            }

            for(const BoundaryLink &bd: cell.boundaries())
            {
                const Vector2D& rf = bd.rFaceVec();
                const Vector2D& sf = bd.outwardNorm();

                grad = outer(field(bd.face()) - field(cell), rf/dot(rf, rf));

                gradField(cell) += Tensor2D(grad.xx*fabs(sf.x), grad.xy*fabs(sf.y),
                                            grad.yx*fabs(sf.x), grad.yy*fabs(sf.y));

                sumSfx += fabs(sf.x);
                sumSfy += fabs(sf.y);
            }

            gradField(cell) = Tensor2D(gradField(cell).xx/sumSfx, gradField(cell).xy/sumSfy,
                                       gradField(cell).yx/sumSfx, gradField(cell).yy/sumSfy);
        }

        for(const Face& face: field.grid.interiorFaces())
        {
            const Vector2D rc = face.rCell().centroid() - face.lCell().centroid();
            gradField(face) = outer(field(face.rCell()) - field(face.lCell()), rc/dot(rc, rc));
        }

        for(const Face& face: field.grid.boundaryFaces())
        {
            const Vector2D rf = face.centroid() - face.lCell().centroid();
            gradField(face) = outer(field(face) - field(face.lCell()), rf/dot(rf, rf));
        }

        break;

    default:
        throw Exception("fv", "computeGradient", "unrecognized gradient evaluation method.");
    }
}

}
