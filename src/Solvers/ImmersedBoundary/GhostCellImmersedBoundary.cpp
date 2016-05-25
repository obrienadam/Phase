#include "GhostCellImmersedBoundary.h"
#include "BilinearInterpolation.h"

namespace gc {

Equation<ScalarFiniteVolumeField> ib(const ImmersedBoundaryObject& ibObj, ScalarFiniteVolumeField& field)
{
    typedef Equation<ScalarFiniteVolumeField>::Triplet Triplet;

    std::vector<Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);

    const UniqueCellGroup &ibCells = field.grid.cellGroup("ibCells");

    entries.reserve(5*ibCells.size());

    for(const ImmersedBoundaryObject::ImmersedBoundaryStencil &ibStencil: ibObj.stencils())
    {
        const size_t row = ibStencil.cell().globalIndex();

        const std::vector< Ref<const Cell> > &kNN = ibStencil.kNearestNeighbours();
        std::vector<Point2D> centroids = {
            kNN[0].get().centroid(),
            kNN[1].get().centroid(),
            kNN[2].get().centroid(),
            kNN[3].get().centroid(),
        };

        BilinearInterpolation bi(centroids);
        std::vector<Scalar> coeffs = bi(ibStencil.imagePoint());

        int cols[] = {
            kNN[0].get().globalIndex(),
            kNN[1].get().globalIndex(),
            kNN[2].get().globalIndex(),
            kNN[3].get().globalIndex(),
        };

        Scalar centralCoeff, gammaIp, gamma[4];
        Vector2D nWall, gradGamma, nl;
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

            gamma[0] = ibObj.csf().gamma()[kNN[0].get().id()];
            gamma[1] = ibObj.csf().gamma()[kNN[1].get().id()];
            gamma[2] = ibObj.csf().gamma()[kNN[2].get().id()];
            gamma[3] = ibObj.csf().gamma()[kNN[3].get().id()];

            gammaIp = bi(gamma, ibStencil.imagePoint());
            nWall = ibStencil.cell().centroid() - ibStencil.imagePoint();

            gradGamma = (ibObj.csf().gamma()[ibStencil.cell().id()] - gammaIp)*nWall.unitVec();

            nl = ibObj.csf().computeContactLineNormal(gradGamma, nWall);

            gradGamma = -gradGamma.mag()*nl;

            eqn.boundaries()(row) -= nWall.mag()*dot(gradGamma, nWall.unitVec());
            centralCoeff = -1.;

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

    const UniqueCellGroup &ibCells = field.grid.cellGroup("ibCells");

    entries.reserve(10*ibCells.size());

    for(const ImmersedBoundaryObject::ImmersedBoundaryStencil &ibStencil: ibObj.stencils())
    {
        const size_t rowX = ibStencil.cell().globalIndex();
        const size_t rowY = rowX + nActiveCells;

        const std::vector< Ref<const Cell> > &kNN = ibStencil.kNearestNeighbours();
        std::vector<Point2D> centroids = {
            kNN[0].get().centroid(),
            kNN[1].get().centroid(),
            kNN[2].get().centroid(),
            kNN[3].get().centroid(),
        };

        BilinearInterpolation bi(centroids);
        std::vector<Scalar> coeffs = bi(ibStencil.imagePoint());

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
        Scalar centralCoeff;
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
