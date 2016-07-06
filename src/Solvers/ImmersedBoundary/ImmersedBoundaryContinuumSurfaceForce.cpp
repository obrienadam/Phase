#include "ImmersedBoundaryContinuumSurfaceForce.h"
#include "ImmersedBoundaryObject.h"
#include "BilinearInterpolation.h"

ImmersedBoundaryContinuumSurfaceForce::ImmersedBoundaryContinuumSurfaceForce(const Input &input,
                                                                             const ScalarFiniteVolumeField &gamma,
                                                                             const VectorFiniteVolumeField &u,
                                                                             std::map<std::string, ScalarFiniteVolumeField> &scalarFields,
                                                                             std::map<std::string, VectorFiniteVolumeField> &vectorFields)
    :
      ContinuumSurfaceForce(input, gamma, u, scalarFields, vectorFields)
{

}

VectorFiniteVolumeField ImmersedBoundaryContinuumSurfaceForce::compute(const std::vector<ImmersedBoundaryObject> &ibObjs)
{
    computeGradGammaTilde();
    computeInterfaceNormals(ibObjs);
    computeCurvature();

    VectorFiniteVolumeField ft(gamma_.grid, "ft");
    VectorFiniteVolumeField gradGamma = grad(gamma_);

    for(const Cell &cell: gamma_.grid.fluidCells())
        ft[cell.id()] = sigma_*kappa_[cell.id()]*gradGamma[cell.id()];

    return ft;
}

//- Protected methods

void ImmersedBoundaryContinuumSurfaceForce::computeInterfaceNormals(const std::vector<ImmersedBoundaryObject>& ibObjs)
{
    VectorEquation eqn(n_, "IB contact line normal", SparseMatrix::NoPreconditioner);
    const Scalar centralCoeff = 1.;

    for(const Cell &cell: n_.grid.fluidCells())
    {
        size_t rowX = cell.globalIndex();
        size_t rowY = rowX + n_.grid.nActiveCells();

        eqn.matrix().insert(rowX, rowX) = 1.;
        eqn.matrix().insert(rowY, rowY) = 1.;

        Vector2D n = gradGammaTilde_[cell.id()] == Vector2D(0., 0.) ? Vector2D(0., 0.) : -gradGammaTilde_[cell.id()].unitVec();

        eqn.sources()(rowX) = n.x;
        eqn.sources()(rowY) = n.y;
    }

    for(const ImmersedBoundaryObject &ibObj: ibObjs)
        for(const Cell &cell: ibObj.cells())
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

            //- From previous values, workout the direction and orientation of the contact line
            std::vector<Vector2D> velVals;
            for(const Cell &cell: kNN)
                velVals.push_back(u_[cell.id()]);

            std::vector<Vector2D> gradGammaVals;
            for(const Cell &cell: kNN)
                gradGammaVals.push_back(gradGammaTilde_[cell.id()]);

            const Vector2D uIp = bi(velVals, imagePoint);
            const Vector2D gradGammaIp = bi(gradGammaVals, imagePoint);

            Vector2D bpNormal = gradGammaIp == Vector2D(0., 0.) ? Vector2D(0., 0.) : SurfaceTensionForce::computeContactLineNormal(gradGammaIp, cell.centroid() - imagePoint, uIp);

            eqn.matrix().insert(rowX, rowX) = centralCoeff/2.;
            eqn.matrix().insert(rowY, rowY) = centralCoeff/2.;

            for(int i = 0; i < coeffs.size(); ++i)
            {
                eqn.matrix().coeffRef(rowX, cols[i]) = coeffs[i]/2.;
                eqn.matrix().coeffRef(rowY, cols[i] + n_.grid.nActiveCells()) = coeffs[i]/2.;
            }

            eqn.sources()(rowX) = bpNormal.x;
            eqn.sources()(rowY) = bpNormal.y;
        }

    Scalar error = eqn.solve();

    if(isnan(error))
    {
        printf("Warning: failed to solve the IB normal equation. Attempting to compute normals using ContinuumSurfaceForce::computerInterfaceNormals.\n");
        ContinuumSurfaceForce::computeInterfaceNormals();
        return;
    }

    interpolateFaces(n_);

    for(const Cell &cell: n_.grid.activeCells())
        n_[cell.id()] = n_[cell.id()] == Vector2D(0., 0.) ? Vector2D(0., 0.) : n_[cell.id()].unitVec();

    for(Vector2D &n: n_.faces())
        n = n == Vector2D(0., 0.) ? Vector2D(0., 0.) : n.unitVec();
}
