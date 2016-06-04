#include "ImmersedBoundaryContinuumSurfaceForce.h"
#include "ImmersedBoundaryObject.h"
#include "BilinearInterpolation.h"

ImmersedBoundaryContinuumSurfaceForce::ImmersedBoundaryContinuumSurfaceForce(const Input &input, const ScalarFiniteVolumeField &gamma, const ScalarFiniteVolumeField &rho)
    :
      ContinuumSurfaceForce(input, gamma, rho)
{

}

VectorFiniteVolumeField ImmersedBoundaryContinuumSurfaceForce::compute(const ImmersedBoundaryObject &ibObj)
{
    computeGradGammaTilde();
    computeInterfaceNormals(ibObj);
    // ContinuumSurfaceForce::computeInterfaceNormals();
    computeCurvature();

    VectorFiniteVolumeField ft(gamma_.grid, "ft");

    for(const Cell &cell: gamma_.grid.fluidCells())
        ft[cell.id()] = sigma_*kappa_[cell.id()]*gradGammaTilde_[cell.id()]*rho_[cell.id()]/avgRho_;

    return ft;
}

//- Protected methods

void ImmersedBoundaryContinuumSurfaceForce::computeInterfaceNormals(const ImmersedBoundaryObject& ibObj)
{
    VectorEquation eqn(n_, "IB contact line normal", SparseMatrix::IncompleteLUT);
    const Scalar centralCoeff = 1.;

    for(const Cell &cell: n_.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + n_.grid.nActiveCells();

        eqn.matrix().insert(rowX, rowX) = 1.;
        eqn.matrix().insert(rowY, rowY) = 1.;

        Vector2D n = gradGammaTilde_[cell.id()] == Vector2D(0., 0.) ? Vector2D(0., 0.) : gradGammaTilde_[cell.id()].unitVec();

        eqn.sources()(rowX) = n.x;
        eqn.sources()(rowY) = n.y;
    }

    for(const Cell &cell: n_.grid.cellGroup("ibCells"))
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + n_.grid.nActiveCells();

        Point2D imagePoint = ibObj.imagePoint(cell.centroid());
        std::vector< Ref<const Cell> > kNN = ibObj.boundingCells(imagePoint);

        std::vector<Point2D> centroids = {
            kNN[0].get().centroid(),
            kNN[1].get().centroid(),
            kNN[2].get().centroid(),
            kNN[3].get().centroid(),
        };

        std::vector<int> cols = {
            kNN[0].get().globalIndex(),
            kNN[1].get().globalIndex(),
            kNN[2].get().globalIndex(),
            kNN[3].get().globalIndex(),
        };

        BilinearInterpolation bi(centroids);
        std::vector<Scalar> coeffs = bi(imagePoint);

        Vector2D bpNormal = (imagePoint - ibObj.centroid()).x < 0. ? (imagePoint - cell.centroid()).rotate(M_PI - thetaAdv_).unitVec() : (imagePoint - cell.centroid()).rotate(-(M_PI - thetaAdv_)).unitVec();

        eqn.matrix().insert(rowX, rowX) = centralCoeff/2.;
        eqn.matrix().insert(rowY, rowY) = centralCoeff/2.;

        for(int i = 0; i < coeffs.size(); ++i)
        {
            eqn.matrix().insert(rowX, cols[i]) = coeffs[i]/2.;
            eqn.matrix().insert(rowY, cols[i] + n_.grid.nActiveCells()) = coeffs[i]/2.;
        }

        eqn.sources()(rowX) = bpNormal.x;
        eqn.sources()(rowY) = bpNormal.y;
    }

    eqn.solve();

    for(const Cell &cell: n_.grid.activeCells())
        n_[cell.id()] = n_[cell.id()] == Vector2D(0., 0.) ? Vector2D(0., 0.) : n_[cell.id()].unitVec();

    interpolateFaces(n_);

    for(Vector2D &n: n_.faces())
        n = n == Vector2D(0., 0.) ? Vector2D(0., 0.) : n.unitVec();
}
