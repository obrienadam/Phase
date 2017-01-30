#include "Poisson.h"
#include "Laplacian.h"
#include "EigenSparseMatrixSolver.h"

Poisson::Poisson(const Input &input, const Communicator &comm, FiniteVolumeGrid2D &grid)
    :
      Solver(input, comm, grid),
      phi(addScalarField(input, "phi")),
      gamma(addScalarField("gamma")),
      phiEqn_(input, comm, phi, "phiEqn")
{
    //- All active cells to fluid cells
    grid_.createCellZone("fluid", grid_.getCellIds(grid_.activeCells()));
    ib_.initCellZones();
    grid_.computeOrdering(comm);

    gamma.fill(input.caseInput().get<Scalar>("Properties.gamma", 1.));
    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);
}

Scalar Poisson::solve(Scalar timeStep)
{
    phiEqn_ = (fv::laplacian(gamma, phi) + ib_.eqns(phi) == 0.);
    Scalar error = phiEqn_.solve();

    grid_.sendMessages(comm_, phi);

    return error;
}
