#include "ContinuumSurfaceForce.h"

ContinuumSurfaceForce::ContinuumSurfaceForce(const Input &input, const ScalarFiniteVolumeField &gamma)
    :
      SurfaceTensionForce(input, gamma),
      gammaTilde_(gamma.grid, "gammaTilde"),
      n_(gamma.grid, "interfaceNormals"),
      kappa_(gamma.grid, "kappa")
{
    CellSearch cs(gamma_.grid.activeCells());
    kernelWidth_ = input.caseInput().get<Scalar>("Solver.smoothingKernelRadius");

    cellRangeSearch_.clear();
    cellRangeSearch_.resize(gamma_.grid.cells().size());

    for(const Cell &cell: gamma_.grid.cells())
        cellRangeSearch_[cell.id()] = cs.rangeSearch(Circle(cell.centroid(), kernelWidth_));
}

VectorFiniteVolumeField ContinuumSurfaceForce::compute()
{
    computeInterfaceNormals();
    computeCurvature();

    return VectorFiniteVolumeField(sigma_*kappa_*grad(gammaTilde_));
}

//- Private methods

void ContinuumSurfaceForce::computeInterfaceNormals()
{
    gammaTilde_ = smooth(gamma_, cellRangeSearch_, kernelWidth_);
    interpolateFaces(gammaTilde_);
    n_ = grad(gammaTilde_);

    for(Vector2D &vec: n_)
    {
        if(vec == Vector2D(0., 0.))
            vec.x = vec.y = 0.;
        else
            vec = vec.unitVec();

        if(isnan(vec.magSqr()))
            throw Exception("ContinuumSurfaceFoce", "computeInterfaceNormals", "NaN value detected.");
    }

    interpolateFaces(n_);
}

void ContinuumSurfaceForce::computeCurvature()
{
    for(const Cell &cell: kappa_.grid.fluidCells())
    {
        Scalar &k = kappa_[cell.id()] = 0.;

        for(const InteriorLink &nb: cell.neighbours())
            k -= dot(n_.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
        {
            k -= dot(n_.faces()[bd.face().id()], bd.outwardNorm()); // Modify surface tension force near boundary
        }

        k /= cell.volume();

        if(isnan(k))
            k = 0.;
            //throw Exception("ContinuumSurfaceFoce", "computeCurvature", "NaN value detected. Cell volume = " + std::to_string(cell.volume()));
    }
}
