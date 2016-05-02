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

        std::vector<Scalar> coeffs = BilinearInterpolation(centroids)(ibStencil.imagePoint());

        int cols[] = {
            kNN[0].get().globalIndex(),
            kNN[1].get().globalIndex(),
            kNN[2].get().globalIndex(),
            kNN[3].get().globalIndex(),
        };

        Scalar centralCoeff;
        // Then here we shall do the boundary assembly!
        switch(ibObj.boundaryType(field.name))
        {
        case ImmersedBoundaryObject::FIXED:
            centralCoeff = 1.;
            break;
        case ImmersedBoundaryObject::NORMAL_GRADIENT:
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
