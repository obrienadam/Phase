#include "FiniteVolume/Discretization/AxisymmetricTimeDerivative.h"
#include "FiniteVolume/Discretization/AxisymmetricDivergence.h"
#include "FiniteVolume/Discretization/AxisymmetricExplicitDivergence.h"
#include "FiniteVolume/Discretization/AxisymmetricLaplacian.h"
#include "FiniteVolume/Discretization/AxisymmetricStressTensor.h"
#include "FiniteVolume/Discretization/AxisymmetricSource.h"

#include "FractionalStepAxisymmetric.h"

FractionalStepAxisymmetric::FractionalStepAxisymmetric(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStep(input, grid)
{

}

Scalar FractionalStepAxisymmetric::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for (const Cell &cell: *fluid_)
    {
        Scalar co = 0.;

        for (const InteriorLink &nb: cell.neighbours())
            co += std::max(dot(u_(nb.face()), nb.polarOutwardNorm()), 0.);

        for (const BoundaryLink &bd: cell.boundaries())
            co += std::max(dot(u_(bd.face()), bd.polarOutwardNorm()), 0.);

        co *= timeStep / cell.polarVolume();
        co_(cell) = co;
        maxCo = std::max(co, maxCo);
    }

    co_.sendMessages();

    return grid_->comm().max(maxCo);
}

Scalar FractionalStepAxisymmetric::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 2);
    uEqn_ = (axi::ddt(u_, timeStep) + axi::dive(u_, u_, 0.5)
             == axi::laplacian(mu_ / rho_, u_, 0.5) - axi::src::src(gradP_));

    Scalar error = uEqn_.solve();

    for(const Cell &c: *fluid_)
        u_(c) += timeStep * gradP_(c);

    u_.sendMessages();
    u_.interpolateFaces();

    return error;
}

Scalar FractionalStepAxisymmetric::solvePEqn(Scalar timeStep)
{
    pEqn_ = (axi::laplacian(timeStep, p_) == axi::src::div(u_));
    Scalar error = pEqn_.solve();
    p_.sendMessages();
    p_.setBoundaryFaces();

    gradP_.computeAxisymmetric(*fluid_);
    gradP_.sendMessages();

    return error;
}

void FractionalStepAxisymmetric::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: grid_->localCells())
        u_(cell) -= timeStep * gradP_(cell);

    for (const Face &face: grid_->faces())
        u_(face) -= timeStep * gradP_(face);

    u_.sendMessages();
}

Scalar FractionalStepAxisymmetric::maxDivergenceError()
{
    Scalar maxError = 0.;

    for (const Cell &cell: *fluid_)
    {
        Scalar divU = 0.;
        for (const InteriorLink &nb: cell.neighbours())
            divU += dot(u_(nb.face()), nb.face().polarOutwardNorm(cell.centroid()));

        for (const BoundaryLink &bd: cell.boundaries())
            divU += dot(u_(bd.face()), bd.face().polarOutwardNorm(cell.centroid()));

        maxError = std::abs(divU) > maxError ? std::abs(divU) : maxError;
    }

    return grid_->comm().max(maxError);
}
