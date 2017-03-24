#include <functional>

#include "FaceInterpolation.h"

namespace fv
{

    template<class T>
    void interpolateFaces(InterpolationMethod method, FiniteVolumeField<T> &field)
    {
        std::function<Scalar(const Face &)> alpha;

        switch (method)
        {
            case UNWEIGHTED:
                alpha = [](const Face &face) { return 0.5; };
                break;

            case INVERSE_VOLUME:
                alpha = [](const Face &face) {
                    return face.rCell().volume() / (face.lCell().volume() + face.rCell().volume());
                };
                break;

            case INVERSE_DISTANCE:
                alpha = [](const Face &face) {
                    const Scalar rd = (face.centroid() - face.rCell().centroid()).mag();
                    const Scalar ld = (face.centroid() - face.lCell().centroid()).mag();

                    return rd / (ld + rd);
                };

                break;

            case INVERSE_SQR_DISTANCE:
                alpha = [](const Face &face) {
                    const Scalar rd = (face.centroid() - face.rCell().centroid()).magSqr();
                    const Scalar ld = (face.centroid() - face.lCell().centroid()).magSqr();

                    return rd / (ld + rd);
                };

                break;

        };

        for (const Face &face: field.grid.interiorFaces())
        {
            const Scalar tmp = alpha(face);

            field(face) = field(face.lCell()) * tmp + field(face.rCell()) * (1. - tmp);
        }

        for (const Face &face: field.grid.boundaryFaces())
        {
            switch (field.boundaryType(face))
            {
                case FiniteVolumeField<T>::FIXED:
                    break;

                default:
                    field(face) = field(face.lCell());
            }
        }
    }

    template<class T>
    void upwindInterpolateFaces(const FiniteVolumeField<Vector2D>& u, FiniteVolumeField<T> &field)
    {
        for(const Face &face: field.grid.interiorFaces())
        {
            Scalar flux = dot(u(face), face.outwardNorm(face.lCell().centroid()));
            field(face) = flux > 0. ? field(face.lCell()) : field(face.rCell());
        }

        for (const Face &face: field.grid.boundaryFaces())
        {
            switch (field.boundaryType(face))
            {
                case FiniteVolumeField<T>::FIXED:
                    break;

                default:
                    field(face) = field(face.lCell());
            }
        }
    }

    template<class T>
    void harmonicInterpolateFaces(InterpolationMethod method, FiniteVolumeField<T> &field)
    {
        std::function<Scalar(const Face &)> alpha;

        switch (method)
        {
            case UNWEIGHTED:
                alpha = [](const Face &face) { return 0.5; };
                break;

            case INVERSE_VOLUME:
                alpha = [](const Face &face) {
                    return face.rCell().volume() / (face.lCell().volume() + face.rCell().volume());
                };
                break;

            case INVERSE_DISTANCE:
                alpha = [](const Face &face) {
                    const Scalar rd = (face.centroid() - face.rCell().centroid()).mag();
                    const Scalar ld = (face.centroid() - face.lCell().centroid()).mag();

                    return rd / (ld + rd);
                };

                break;

            case INVERSE_SQR_DISTANCE:
                alpha = [](const Face &face) {
                    const Scalar rd = (face.centroid() - face.rCell().centroid()).magSqr();
                    const Scalar ld = (face.centroid() - face.lCell().centroid()).magSqr();

                    return rd / (ld + rd);
                };

                break;

        };

        for (const Face &face: field.grid.interiorFaces())
        {
            const Scalar tmp = alpha(face);
            field(face) = 1. / (tmp / field(face.lCell()) + (1. - tmp) / field(face.rCell()));
        }

        for (const Face &face: field.grid.boundaryFaces())
        {
            switch (field.boundaryType(face))
            {
                case FiniteVolumeField<T>::FIXED:
                    break;

                default:
                    field(face) = field(face.lCell());
            }
        }
    }

    template<class T>
    void interpolateFaces(const FiniteVolumeField<Scalar> &w, FiniteVolumeField<T> &field, bool inverseWeighting)
    {
        if (inverseWeighting)
        {
            for (const Face &face: field.grid.interiorFaces())
                field(face) = (field(face.lCell()) / w(face.lCell()) + field(face.rCell()) / w(face.rCell())) /
                              (1. / w(face.lCell()) + 1. / w(face.rCell()));
        }
        else
        {
            for (const Face &face: field.grid.interiorFaces())
                field(face) = (w(face.lCell()) * field(face.lCell()) + w(face.rCell()) * field(face.rCell())) /
                              (w(face.lCell()) + w(face.rCell()));
        }

        for (const Face &face: field.grid.boundaryFaces())
        {
            switch (field.boundaryType(face))
            {
                case FiniteVolumeField<T>::FIXED:
                    break;

                default:
                    field(face) = field(face.lCell());
            }
        }
    }

}
