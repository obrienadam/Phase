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
            centralCoeff = 1;
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

        if(ibObj.boundaryType(field.name) == ImmersedBoundaryObject::CONTACT_ANGLE)
        {
            Scalar theta = ibObj.csf().theta();
            const Point2D boundaryPoint = ibObj.boundaryPoint(cell.centroid());

            Scalar rotAngle = (boundaryPoint - ibObj.centroid()).x > 0. ? theta - M_PI/2. : M_PI/2. - theta;

            Vector2D ip2 = (imagePoint - boundaryPoint).rotate(rotAngle) + imagePoint;

            kNN = ibObj.boundingCells(ip2);

            centroids = {
                        kNN[0].get().centroid(),
                        kNN[1].get().centroid(),
                        kNN[2].get().centroid(),
                        kNN[3].get().centroid(),
                    };

            bi = BilinearInterpolation(centroids);

            std::vector<Scalar> tmpCoeffs = bi(ip2);

            for(Scalar coeff: tmpCoeffs)
                coeffs.push_back(-2*coeff);

            for(const Cell &cell: kNN)
                cols.push_back(cell.globalIndex());
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
        std::vector<Scalar> coeffs = bi(imagePoint);

        int colsX[] = {
            kNN[0].get().globalIndex(),
            kNN[1].get().globalIndex(),
            kNN[2].get().globalIndex(),
            kNN[3].get().globalIndex(),
        };

        int colsY[4];
        for(int i = 0; i < 4; ++i)
            colsY[i] = colsX[i] + nActiveCells;

        // Then here we shall do the boundary assembly!
        Scalar centralCoeff, l;
        switch(ibObj.boundaryType(field.name))
        {
        case ImmersedBoundaryObject::FIXED:
            centralCoeff = 1.;
            break;
        case ImmersedBoundaryObject::NORMAL_GRADIENT:
            centralCoeff = -1.;
            break;

        case ImmersedBoundaryObject::CONTACT_ANGLE:
            throw Exception("gc", "ib", "contact angle boundaries are not valid for vector fields.");

        case ImmersedBoundaryObject::PARTIAL_SLIP:
            l = (cell.centroid() - imagePoint).mag();
            centralCoeff = 1. + 2.*ibObj.boundaryRefValue(field.name)/l;

            for(Scalar &coeff: coeffs)
                coeff *= 1. - 2.*ibObj.boundaryRefValue(field.name)/l;

            break;

        default:
            throw Exception("gc", "ib", "invalid boundary type.");
        }

        entries.push_back(Triplet(rowX, rowX, centralCoeff));
        entries.push_back(Triplet(rowY, rowY, centralCoeff));

        for(int i = 0; i < coeffs.size(); ++i)
        {
            entries.push_back(Triplet(rowX, colsX[i], coeffs[i]));
            entries.push_back(Triplet(rowY, colsY[i], coeffs[i]));
        }
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

}
