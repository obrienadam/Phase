#include "ContinuumSurfaceForce.h"
#include "BilinearInterpolation.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

ContinuumSurfaceForce::ContinuumSurfaceForce(const Input &input,
                                             Solver &solver)
    :
      SurfaceTensionForce(input, solver),
      gammaTilde_(solver.addScalarField("gammaTilde")),
      gradGammaTilde_(solver.addVectorField("gradGammaTilde"))
{
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");
    curvatureCutoffTolerance_ = input.caseInput().get<Scalar>("Solver.curvatureCutoffTolerance", 1e-14);

    constructSmoothingKernels();
}

VectorFiniteVolumeField ContinuumSurfaceForce::compute()
{
    computeGradient(fv::FACE_TO_CELL, gamma_, gradGamma_);
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    VectorFiniteVolumeField ft(gamma_.grid, "ft");

    for(const Cell &cell: gamma_.grid.fluidCells())
        ft(cell) = sigma_*kappa_(cell)*gradGamma_(cell);

    for(const Face &face: gamma_.grid.faces())
        ft(face) = sigma_*kappa_(face)*gradGamma_(face);

    return ft;
}

void ContinuumSurfaceForce::constructSmoothingKernels()
{
    cellRangeSearch_.clear();
    cellRangeSearch_.resize(gamma_.grid.cells().size());

    for(const Cell &cell: gamma_.grid.activeCells())
        cellRangeSearch_[cell.id()] = gamma_.grid.activeCells().rangeSearch(Circle(cell.centroid(), kernelWidth_));
}

//- Private methods

void ContinuumSurfaceForce::computeGradGammaTilde()
{
    gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    computeGradient(fv::FACE_TO_CELL, gammaTilde_, gradGammaTilde_);
}

void ContinuumSurfaceForce::computeInterfaceNormals()
{
    VectorEquation eqn(n_, "IB contact line normal");
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

    for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
        for(const Cell &cell: ibObj.cells())
        {
            size_t rowX = cell.globalIndex();
            size_t rowY = rowX + n_.grid.nActiveCells();

            const std::vector< Ref<const Cell> > &kNN = ibObj.imagePointCells(cell);
            const BilinearInterpolation &bi = ibObj.imagePointInterpolation(cell);
            const Point2D &imagePoint = ibObj.imagePoint(cell);

            std::vector<Scalar> coeffs = bi(imagePoint);

            std::vector<Index> cols = {
                kNN[0].get().globalIndex(),
                kNN[1].get().globalIndex(),
                kNN[2].get().globalIndex(),
                kNN[3].get().globalIndex(),
            };

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

    interpolateFaces(fv::INVERSE_VOLUME, n_);

    for(const Cell &cell: n_.grid.activeCells())
        n_[cell.id()] = n_[cell.id()] == Vector2D(0., 0.) ? Vector2D(0., 0.) : n_[cell.id()].unitVec();

    for(Vector2D &n: n_.faces())
        n = n == Vector2D(0., 0.) ? Vector2D(0., 0.) : n.unitVec();
}

void ContinuumSurfaceForce::computeCurvature()
{
    for(const Cell &cell: kappa_.grid.fluidCells())
    {
        Scalar &k = kappa_[cell.id()] = 0.;

        for(const InteriorLink &nb: cell.neighbours())
            k += dot(n_.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
            k += dot(n_.faces()[bd.face().id()], bd.outwardNorm());

        k /= cell.volume();
    }

    interpolateFaces(fv::INVERSE_VOLUME, kappa_);
}

void ContinuumSurfaceForce::interpolateCurvatureFaces()
{
    for(const Face &face: kappa_.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();
        const Vector2D rc = rCell.centroid() - lCell.centroid();

        const Scalar gradGammaMagSqr = ((gamma_(rCell) - gamma_(lCell))*rc/dot(rc, rc)).magSqr();

        if(gradGammaMagSqr < curvatureCutoffTolerance_)
            kappa_(face) = 0.;
        else
        {
            const Scalar g = rCell.volume()/(rCell.volume() + lCell.volume());
            kappa_(face) = g*kappa_(lCell) + (1. - g)*kappa_(rCell);
        }
    }

    for(const Face &face: kappa_.grid.boundaryFaces())
    {
        const Cell& cell = face.lCell();
        const Vector2D rf = face.centroid() - cell.centroid();

        const Scalar gradGammaMagSqr = ((gamma_(face) - gamma_(cell))*rf/dot(rf, rf)).magSqr();

        if(gradGammaMagSqr < curvatureCutoffTolerance_)
            kappa_(face) = 0.;
        else
            kappa_(face) = kappa_(cell);
    }
}
