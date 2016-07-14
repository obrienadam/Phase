#include "FractionalStepMultiphase.h"
#include "Cicsam.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"

FractionalStepMultiphase::FractionalStepMultiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      FractionalStep(grid, input),
      gamma(addScalarField(input, "gamma")),
      ft(addVectorField("ft")),
      gammaEqn_(gamma, "gammaEqn"),
      surfaceTensionForce_(input, gamma, u, scalarFields_, vectorFields_)
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    setInitialConditions(input);

    //surfaceTensionForce_->compute();
    computeRho();
    computeMu();
}

Scalar FractionalStepMultiphase::solve(Scalar timeStep)
{
    computeRho();
    computeMu();

    u.savePreviousTimeStep(timeStep, 1);
    gamma.savePreviousTimeStep(timeStep, 1);

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    solveGammaEqn(timeStep);

    return 0.;
}

//- Protected methods

Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    ft = surfaceTensionForce_.compute();

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u) + ib_.eqns(u) == ab::laplacian(mu, u) + fv::source(ft));
    Scalar error = uEqn_.solve();
    interpolateFaces(u);
    return error;
}

Scalar FractionalStepMultiphase::solvePEqn(Scalar timeStep)
{
    divUStar.fill(0.);

    for(const Cell &cell: grid_.fluidCells())
    {
        for(const InteriorLink &nb: cell.neighbours())
            divUStar[cell.id()] += rho[cell.id()]/timeStep*dot(u.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
            divUStar[cell.id()] += rho[cell.id()]/timeStep*dot(u.faces()[bd.face().id()], bd.outwardNorm());
    }

    pEqn_ = (fv::laplacian(p) + ib_.eqns(p) == divUStar);
    Scalar error = pEqn_.solve();

    interpolateFaces(p);

    return error;
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    interpolateFaces(gamma);

    gammaEqn_ = (cicsam::div(u, gamma, timeStep, cicsam::HC) == 0.);
    Scalar error = gammaEqn_.solve();
    return error;
}

void FractionalStepMultiphase::computeRho()
{
    const ScalarFiniteVolumeField &gammaTilde = surfaceTensionForce_.gammaTilde();

    for(const Cell &cell: grid_.activeCells())
        rho[cell.id()] = (1. - gammaTilde[cell.id()])*rho1_ + gammaTilde[cell.id()]*rho2_;

    harmonicInterpolateFaces(rho);
}

void FractionalStepMultiphase::computeMu()
{
    const ScalarFiniteVolumeField &gammaTilde = surfaceTensionForce_.gammaTilde();

    for(const Cell &cell: grid_.activeCells())
        mu[cell.id()] = (1. - gammaTilde[cell.id()])*mu1_ + gammaTilde[cell.id()]*mu2_;

    harmonicInterpolateFaces(mu);
}
