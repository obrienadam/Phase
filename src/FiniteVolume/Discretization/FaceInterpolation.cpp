#include "FaceInterpolation.h"
#include "VectorFiniteVolumeField.h"

namespace fv
{

void interpolateFaces(InterpolationMethod method, VectorFiniteVolumeField& field)
{
    std::function<Scalar(const Face&)> alpha;

    switch(method)
    {
    case UNWEIGHTED:
        alpha = [](const Face& face) { return 0.5; };
        break;

    case INVERSE_VOLUME:
        alpha = [](const Face& face) { return face.rCell().volume()/(face.lCell().volume() + face.rCell().volume()); };
        break;

    case INVERSE_DISTANCE:
        alpha = [](const Face& face)
        {
            const Scalar rd = (face.centroid() - face.rCell().centroid()).mag();
            const Scalar ld = (face.centroid() - face.lCell().centroid()).mag();

            return rd/(ld + rd);
        };

        break;

    case INVERSE_SQR_DISTANCE:
        alpha = [](const Face& face)
        {
            const Scalar rd = (face.centroid() - face.rCell().centroid()).magSqr();
            const Scalar ld = (face.centroid() - face.lCell().centroid()).magSqr();

            return rd/(ld + rd);
        };

        break;

    };

    for(const Face& face: field.grid.interiorFaces())
    {
        const Scalar tmp = alpha(face);

        field(face) = field(face.lCell())*tmp + field(face.rCell())*(1. - tmp);
    }

    for(const Face& face: field.grid.boundaryFaces())
    {
        switch(field.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT: case VectorFiniteVolumeField::OUTFLOW:
            field(face) = field(face.lCell());
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            field(face) = field(face.lCell()) - dot(field(face.lCell()), nWall)*nWall/nWall.magSqr();
            break;
        }

        default:
            throw Exception("VectorFiniteVolumeField", "interpolateFaces", "unrecongnized boundary condition type.");
        }
    }
}

}
