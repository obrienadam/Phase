#include <functional>

#include "FaceInterpolation.h"

namespace fv
{

template<class T>
void interpolateFaces(InterpolationMethod method, FiniteVolumeField<T>& field)
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

        field.faces()[face.id()] = field[face.lCell().id()]*tmp + field[face.rCell().id()]*(1. - tmp);
    }

    for(const Face& face: field.grid.boundaryFaces())
    {
        switch(field.boundaryType(face.id()))
        {
        case FiniteVolumeField<T>::FIXED:
            break;

        case FiniteVolumeField<T>::NORMAL_GRADIENT: case FiniteVolumeField<T>::OUTFLOW:
            field.faces()[face.id()] = field[face.lCell().id()];
            break;

        case FiniteVolumeField<T>::SYMMETRY:
            field.faces()[face.id()] = field[face.lCell().id()];
            break;

        default:
            throw Exception("fv", "interpolateFaces", "unrecongnized boundary condition type.");
        }
    }
}

template<class T>
void harmonicInterpolateFaces(InterpolationMethod method, FiniteVolumeField<T>& field)
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
        field(face) = 1./(tmp/field(face.lCell()) + (1. - tmp)/field(face.rCell()));
    }

    for(const Face& face: field.grid.boundaryFaces())
    {
        switch(field.boundaryType(face.id()))
        {
        case FiniteVolumeField<T>::FIXED:
            break;

        case FiniteVolumeField<T>::NORMAL_GRADIENT:
            field(face) = field(face.lCell());
            break;

        case FiniteVolumeField<T>::SYMMETRY:
            field(face) = field(face.lCell());
            break;

        default:
            throw Exception("fv", "harmonicInterpolateFaces", "unrecongnized boundary condition type.");
        }
    }
}

template<class T>
void interpolateFaces(const FiniteVolumeField<Scalar> &w, FiniteVolumeField<T>& field, bool inverseWeighting)
{
    if(inverseWeighting)
    {
        for(const Face& face: field.grid.interiorFaces())
            field(face) = (field(face.lCell())/w(face.lCell()) + field(face.rCell())/w(face.rCell()))/(1./w(face.lCell()) + 1./w(face.rCell()));
    }
    else
    {
        for(const Face& face: field.grid.interiorFaces())
            field(face) = (w(face.lCell())*field(face.lCell()) + w(face.rCell())*field(face.rCell()))/(w(face.lCell()) + w(face.rCell()));
    }

    for(const Face& face: field.grid.boundaryFaces())
    {
        switch(field.boundaryType(face.id()))
        {
        case FiniteVolumeField<T>::FIXED:
            break;

        case FiniteVolumeField<T>::NORMAL_GRADIENT: case FiniteVolumeField<T>::OUTFLOW: case FiniteVolumeField<T>::SYMMETRY:
            field(face) = field(face.lCell());
            break;

        default:
            throw Exception("fv", "interpolateFaces", "unrecongnized boundary condition type.");
        }
    }
}

}
