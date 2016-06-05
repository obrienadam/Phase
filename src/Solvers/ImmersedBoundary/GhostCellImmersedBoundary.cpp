#include "GhostCellImmersedBoundary.h"
#include "BilinearInterpolation.h"

namespace gc {

Equation<ScalarFiniteVolumeField> ib(const ImmersedBoundaryObject& ibObj, ScalarFiniteVolumeField& field)
{
    typedef Equation<ScalarFiniteVolumeField>::Triplet Triplet;

    std::vector<Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    entries.reserve(5*field.grid.nActiveCells());

    for(const Cell &cell: ibObj.cells())
    {
        const size_t row = cell.globalIndex();
        const Point2D imagePoint = ibObj.imagePoint(cell.centroid());

        std::vector< Ref<const Cell> > kNN = ibObj.boundingCells(imagePoint);
        std::vector<Point2D> centroids = {
            kNN[0].get().centroid(),
            kNN[1].get().centroid(),
            kNN[2].get().centroid(),
            kNN[3].get().centroid(),
        };

        BilinearInterpolation bi(centroids);
        std::vector<Scalar> coeffs = bi(imagePoint);

        std::vector<int> cols = {
            kNN[0].get().globalIndex(),
            kNN[1].get().globalIndex(),
            kNN[2].get().globalIndex(),
            kNN[3].get().globalIndex(),
        };

        Scalar centralCoeff, l;
        // Then here we shall do the boundary assembly!
        switch(ibObj.boundaryType(field.name))
        {
        case ImmersedBoundaryObject::FIXED:
            centralCoeff = 1.;
            break;

        case ImmersedBoundaryObject::NORMAL_GRADIENT:
            centralCoeff = -1.;
            break;

        case ImmersedBoundaryObject::CONTACT_ANGLE:
            centralCoeff = -1;
            break;

        case ImmersedBoundaryObject::PARTIAL_SLIP:
            l = (cell.centroid() - imagePoint).mag();
            centralCoeff = 1. + 2.*ibObj.boundaryRefValue(field.name)/l;

            for(Scalar &coeff: coeffs)
                coeff *= 1. - 2.*ibObj.boundaryRefValue(field.name)/l;

            break;

        default:
            throw Exception("gc", "ib", "invalid boundary type.");
        }

        entries.push_back(Triplet(row, row, centralCoeff));
        for(int i = 0; i < coeffs.size(); ++i)
            entries.push_back(Triplet(row, cols[i], coeffs[i]));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}


Equation<VectorFiniteVolumeField> ib(const ImmersedBoundaryObject &ibObj, VectorFiniteVolumeField &field)
{
    typedef Equation<VectorFiniteVolumeField>::Triplet Triplet;

    std::vector<Triplet> entries;
    Equation<VectorFiniteVolumeField> eqn(field);
    const size_t nActiveCells = field.grid.nActiveCells();

    entries.reserve(10*field.grid.nActiveCells());

    for(const Cell &cell: ibObj.cells())
    {
        const size_t rowX = cell.globalIndex();
        const size_t rowY = rowX + nActiveCells;

        const Point2D imagePoint = ibObj.imagePoint(cell.centroid());

        const std::vector< Ref<const Cell> > kNN = ibObj.boundingCells(imagePoint);
        std::vector<Point2D> centroids = {
            kNN[0].get().centroid(),
            kNN[1].get().centroid(),
            kNN[2].get().centroid(),
            kNN[3].get().centroid(),
        };

        BilinearInterpolation bi(centroids);
        std::vector<Scalar> coeffsX = bi(imagePoint), coeffsY;
        coeffsY = coeffsX;

        std::vector<int> colsX = {
            kNN[0].get().globalIndex(),
            kNN[1].get().globalIndex(),
            kNN[2].get().globalIndex(),
            kNN[3].get().globalIndex(),
        };

        std::vector<int> colsY(4);
        for(int i = 0; i < 4; ++i)
            colsY[i] = colsX[i] + nActiveCells;

        // Then here we shall do the boundary assembly!
        Scalar centralCoeffX, centralCoeffY, l, tmpScalar;
        Vector2D n, t;
        switch(ibObj.boundaryType(field.name))
        {
        case ImmersedBoundaryObject::FIXED:
            centralCoeffX = centralCoeffY = 1.;
            break;
        case ImmersedBoundaryObject::NORMAL_GRADIENT:
            centralCoeffX = centralCoeffY = -1.;
            break;

        case ImmersedBoundaryObject::CONTACT_ANGLE:
            throw Exception("gc", "ib", "contact angle boundaries are not valid for vector fields.");

        case ImmersedBoundaryObject::PARTIAL_SLIP:

            //- Note: the tangent equation replaces the u.x equation and the normal equation replaces the u.y equation

            //- Enforce the tangent condition
            n = (imagePoint - cell.centroid()).unitVec();
            t = n.tangentVec();

            l = (imagePoint - cell.centroid()).mag();

            tmpScalar = ibObj.boundaryRefValue(field.name);
            tmpScalar = (1. + 2.*tmpScalar/l)/(1. - 2.*tmpScalar/l);

            centralCoeffX = tmpScalar*t.x;

            colsX.push_back(rowY);
            coeffsX.push_back(tmpScalar*t.y);

            for(int i = 0; i < 4; ++i)
            {
                colsX.push_back(colsY[i]);
                coeffsX.push_back(coeffsY[i]*t.y);
            }

            //- Enforce the normal condition
            centralCoeffY = n.y;

            colsY.push_back(rowX);
            coeffsY.push_back(n.x);

            for(int i = 0; i < 4; ++i)
            {
                colsY.push_back(colsX[i]);
                coeffsY.push_back(coeffsX[i]*n.x);
            }

            //- Modify the original coefficients
            for(int i = 0; i < 4; ++i)
            {
                coeffsX[i] *= t.x;
                coeffsY[i] *= n.y;
            }

            break;

        default:
            throw Exception("gc", "ib", "invalid boundary type.");
        }

        entries.push_back(Triplet(rowX, rowX, centralCoeffX));
        entries.push_back(Triplet(rowY, rowY, centralCoeffY));

        for(int i = 0; i < coeffsX.size(); ++i)
        {
            entries.push_back(Triplet(rowX, colsX[i], coeffsX[i]));
            entries.push_back(Triplet(rowY, colsY[i], coeffsY[i]));
        }
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

}
