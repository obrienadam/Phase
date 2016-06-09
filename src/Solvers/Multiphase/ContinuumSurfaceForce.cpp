#include "ContinuumSurfaceForce.h"

ContinuumSurfaceForce::ContinuumSurfaceForce(const Input &input, const ScalarFiniteVolumeField &gamma, const ScalarFiniteVolumeField &rho, std::map<std::string, ScalarFiniteVolumeField> &fields)
    :
      SurfaceTensionForce(input, gamma),
      gradGammaTilde_(gamma.grid, "gammaTilde"),
      rho_(rho),
      gammaTilde_((fields.insert(std::make_pair(std::string("gammaTilde"), ScalarFiniteVolumeField(gamma.grid, "gammaTilde"))).first)->second)
{
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");
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
    interpolateFaces(gammaTilde_);
    gradGammaTilde_ = grad(gammaTilde_);
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
