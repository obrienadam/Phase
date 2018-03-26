#include "Poisson.h"
#include "Laplacian.h"

Poisson::Poisson(const Input &input)
        :
        Solver(input),
        phi(*addScalarField(input, "phi")),
        gamma(*addScalarField("gamma")),
        phiEqn_(input, phi, "phiEqn")
{
    //- All active cells to fluid cells
    grid_->createCellZone("fluid");
    grid_->cellZone("fluid").add(grid_->localActiveCells());

    //- Create ib zones if any
    ib_->initCellZones(grid_->cellZone("solid"));

    gamma.fill(input.caseInput().get<Scalar>("Properties.gamma", 1.));
}

Scalar Poisson::solve(Scalar timeStep)
{
    phiEqn_ = (fv::laplacian(gamma, phi) + ib_->bcs(phi) == 0.);
    Scalar error = phiEqn_.solve();

    grid_->sendMessages(phi);

    return error;
}
