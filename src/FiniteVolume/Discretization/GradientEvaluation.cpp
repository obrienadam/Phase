#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

namespace fv
{

VectorFiniteVolumeField computeGradient(GradientEvaluationMethod method,  ScalarFiniteVolumeField& field, bool useCurrentFacesValues)
{
    VectorFiniteVolumeField gradField(field.grid, "grad_" + field.name);

    computeGradient(method, field, gradField, useCurrentFacesValues);

    return gradField;
}

void computeGradient(GradientEvaluationMethod method, ScalarFiniteVolumeField& field, VectorFiniteVolumeField& gradField, bool useCurrentFaceValues)
{
    gradField.fill(Vector2D(0., 0.));

    switch(method)
    {
    case GREEN_GAUSS_CELL_CENTERED:
        if(!useCurrentFaceValues)
            interpolateFaces(INVERSE_VOLUME, field);

        for(const Cell& cell: field.grid.fluidCells())
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

        for(const Cell& cell: field.grid.fluidCells())
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

        for(const Cell &cell: field.grid.fluidCells())
        {
            Scalar sumSfx = 0., sumSfy = 0.;
            Vector2D grad(0., 0.);
            gradField(cell) = grad;

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

        break;

    default:
        throw Exception("fv", "computeGradient", "unrecognized gradient evaluation method.");
    }
}

}
