#include "FiniteVolume/Discretization/AxisymmetricTimeDerivative.h"
#include "FiniteVolume/Discretization/AxisymmetricDivergence.h"
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
            co += std::max(dot(u_(nb.face()), nb.face().polarOutwardNorm(cell.centroid())), 0.);

        for (const BoundaryLink &bd: cell.boundaries())
            co += std::max(dot(u_(bd.face()), bd.face().polarOutwardNorm(cell.centroid())), 0.);

        co *= timeStep / cell.polarVolume();
        maxCo = std::max(co, maxCo);
    }

    return grid_->comm().max(maxCo);
}

Scalar FractionalStepAxisymmetric::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (axi::ddt(u_, timeStep) + axi::div(u_, u_, 0.)
             == axi::divSigma(mu_ / rho_, p_, u_, 0.5));

    Scalar error = uEqn_.solve();

    for(const Cell &c: *fluid_)
        u_(c) += timeStep * gradP_(c);

    grid_->sendMessages(u_);
    u_.interpolateFaces();

    return error;
}

Scalar FractionalStepAxisymmetric::solvePEqn(Scalar timeStep)
{
    pEqn_ = (axi::laplacian(timeStep, p_) == axi::src::div(u_));
    Scalar error = pEqn_.solve();

    grid_->sendMessages(p_);
    p_.interpolateFaces(ScalarFiniteVolumeField::DISTANCE);
    gradP_.computeAxisymmetric(*fluid_);

    return error;
}

void FractionalStepAxisymmetric::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: grid_->localCells())
        u_(cell) -= timeStep * gradP_(cell);

    grid_->sendMessages(u_);

    for (const Face &face: grid_->faces())
        u_(face) -= timeStep * gradP_(face);
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
