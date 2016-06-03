#include "ContinuumSurfaceForce.h"

ContinuumSurfaceForce::ContinuumSurfaceForce(const Input &input, const ScalarFiniteVolumeField &gamma, const ScalarFiniteVolumeField &rho)
    :
      SurfaceTensionForce(input, gamma),
      gradGammaTilde_(gamma.grid, "gammaTilde"),
      rho_(rho)
{
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");

    cellRangeSearch_.clear();
    cellRangeSearch_.resize(gamma_.grid.cells().size());

    for(const Cell &cell: gamma_.grid.cells())
        cellRangeSearch_[cell.id()] = gamma.grid.activeCells().rangeSearch(Circle(cell.centroid(), kernelWidth_));

    avgRho_ = (input.caseInput().get<Scalar>("Properties.rho1") + input.caseInput().get<Scalar>("Properties.rho2"))/2.;
}

VectorFiniteVolumeField ContinuumSurfaceForce::compute()
{
    computeGradGammaTilde();
    computeInterfaceNormals();
    computeCurvature();

    VectorFiniteVolumeField ft(gamma_.grid, "ft");

    for(const Cell &cell: gamma_.grid.fluidCells())
        ft[cell.id()] = sigma_*kappa_[cell.id()]*gradGammaTilde_[cell.id()]*rho_[cell.id()]/avgRho_;

    return ft;
}

//- Private methods

void ContinuumSurfaceForce::computeGradGammaTilde()
{
    ScalarFiniteVolumeField gammaTilde = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    interpolateFaces(gammaTilde);
    gradGammaTilde_ = grad(gammaTilde);
}

void ContinuumSurfaceForce::computeInterfaceNormals()
{
    n_ = gradGammaTilde_;

    for(Vector2D &vec: n_)
        vec = vec == Vector2D(0., 0.) ? Vector2D(0., 0.) : vec.unitVec();

    interpolateFaces(n_);

    //- Set the contact angle normals
    for(const Patch &patch: contactAnglePatches_)
        for(const Face &face: patch.faces())
        {
            const Cell &cell = face.lCell();
            Vector2D sf = face.outwardNorm(cell.centroid());

            n_.faces()[face.id()] = -computeContactLineNormal(gradGammaTilde_[cell.id()], sf);
        }

    //- Make sure all interpolated faces have the correct magnitude
    for(Vector2D &n: n_.faces())
        n = n == Vector2D(0., 0.) ? Vector2D(0., 0.) : n.unitVec();
}

void ContinuumSurfaceForce::computeCurvature()
{
    for(const Cell &cell: kappa_.grid.fluidCells())
    {
        Scalar &k = kappa_[cell.id()] = 0.;

        for(const InteriorLink &nb: cell.neighbours())
            k -= dot(n_.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
            k -= dot(n_.faces()[bd.face().id()], bd.outwardNorm());

        k /= cell.volume();

        if(isnan(k))
            k = 0.;
    }
}
