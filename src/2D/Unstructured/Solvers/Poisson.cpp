#include "FiniteVolume/Equation/Laplacian.h"

#include "Poisson.h"

Poisson::Poisson(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
        :
        Solver(input, grid),
        solid_(*cells_),
        phi(*addField<Scalar>(input, "phi")),
        phiEqn_(input, phi, "phiEqn")
{
    //- All active cells to solid group
    solid_.add(grid_->localCells());

    //- Create ib zones if any
    //ib_->initCellZones(solid_);

    gamma_ = input.caseInput().get<Scalar>("Properties.gamma", 1.);
}

void Poisson::initialize()
{
    phi.savePreviousTimeStep(0, 1);
}

Scalar Poisson::solve(Scalar timeStep)
{
    phiEqn_ = (fv::laplacian(gamma_, phi, 1) + ib_->bcs(phi) == 0.);
    Scalar error = phiEqn_.solve();

    grid_->sendMessages(phi);

    return error;
}
