#include "ContinuumSurfaceForce.h"
#include "BilinearInterpolation.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"
#include "EigenSparseMatrixSolver.h"

ContinuumSurfaceForce::ContinuumSurfaceForce(const Input &input,
                                             Solver &solver)
    :
      SurfaceTensionForce(input, solver),
      gammaTilde_(solver.addScalarField("gammaTilde")),
      gradGammaTilde_(solver.addVectorField("gradGammaTilde"))
{
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");
    curvatureCutoffTolerance_ = input.caseInput().get<Scalar>("Solver.curvatureCutoffTolerance", 1e-10);

    constructSmoothingKernels();
}

VectorFiniteVolumeField ContinuumSurfaceForce::compute()
{
    computeGradient(fv::FACE_TO_CELL, gamma_, gradGamma_);
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    VectorFiniteVolumeField ft(gamma_.grid, "ft");

    for(const Cell &cell: gamma_.grid.cellZone("fluid"))
        ft(cell) = sigma_*kappa_(cell)*gradGamma_(cell);

    for(const Face &face: gamma_.grid.faces())
        ft(face) = sigma_*kappa_(face)*gradGamma_(face);

    return ft;
}

void ContinuumSurfaceForce::constructSmoothingKernels()
{
    cellRangeSearch_ = gamma_.grid.constructSmoothingKernels(kernelWidth_);
}

//- Private methods

void ContinuumSurfaceForce::computeGradGammaTilde()
{
    gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    computeGradient(fv::FACE_TO_CELL, gammaTilde_, gradGammaTilde_);
}

void ContinuumSurfaceForce::computeInterfaceNormals()
{
    Equation<Vector2D> eqn(n_, "IB contact line normal");
    const Scalar centralCoeff = 1.;

    throw Exception("ContinuumSurfaceForce", "computeInterfaceNormals", "this method has been deprecated and should not be called.");

    for(const Cell &cell: n_.grid.cellZone("fluid"))
    {
        eqn.set(cell, cell, 1.);

        Vector2D n = gradGammaTilde_(cell) == Vector2D(0., 0.) ? Vector2D(0., 0.) : -gradGammaTilde_(cell).unitVec();
        eqn.setSource(cell, n);
    }

    for(const ImmersedBoundaryObject &ibObj: solver_.ib().ibObjs())
        for(const Cell &cell: ibObj.cells())
        {
            const std::vector< Ref<const Cell> > &kNN = ibObj.imagePointCells(cell);
            const BilinearInterpolation &bi = ibObj.imagePointInterpolation(cell);
            const Point2D &imagePoint = ibObj.imagePoint(cell);

            std::vector<Scalar> coeffs = bi(imagePoint);

            //- From previous values, workout the direction and orientation of the contact line
            std::vector<Vector2D> velVals;
            for(const Cell &cell: kNN)
                velVals.push_back(u_(cell));

            std::vector<Vector2D> gradGammaVals;
            for(const Cell &cell: kNN)
                gradGammaVals.push_back(gradGammaTilde_(cell));

            const Vector2D uIp = bi(velVals, imagePoint);
            const Vector2D gradGammaIp = bi(gradGammaVals, imagePoint);

            Vector2D bpNormal = gradGammaIp == Vector2D(0., 0.) ? Vector2D(0., 0.) : SurfaceTensionForce::computeContactLineNormal(gradGammaIp, cell.centroid() - imagePoint, uIp);
            eqn.set(cell, cell, centralCoeff/2.);

            int i = 0;
            for(const Cell& nbCell: kNN)
                eqn.set(cell, nbCell, coeffs[i++]/2.);

            eqn.setSource(cell, bpNormal);
        }

    //Scalar error = eqn.solve(*solver_.newSparseMatrixSolver());
    Scalar error = 0;

    if(std::isnan(error))
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
    for(const Cell &cell: kappa_.grid.cellZone("fluid"))
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
