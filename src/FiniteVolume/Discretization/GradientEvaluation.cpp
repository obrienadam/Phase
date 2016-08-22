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

    default:
        throw Exception("fv", "computeGradient", "unrecognized gradient evaluation method.");
    }
}

}
