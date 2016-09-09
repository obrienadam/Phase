#include "FractionalStepMultiphase.h"
#include "Cicsam.h"
#include "CrankNicolson.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

FractionalStepMultiphase::FractionalStepMultiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      FractionalStep(grid, input),
      gamma(addScalarField(input, "gamma")),
      gradGamma(addVectorField("gradGamma")),
      ft(addVectorField("ft")),
      gammaEqn_(input, gamma, "gammaEqn"),
      surfaceTensionForce_(input, *this)
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);

    setInitialConditions(input);

    surfaceTensionForce_.compute();
    computeRho();
    computeMu();
}

Scalar FractionalStepMultiphase::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    solveGammaEqn(timeStep);

    printf("Max Co = %lf\n", courantNumber(timeStep));

    return 0.;
}

//- Protected methods

Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    sg.savePreviousTimeStep(timeStep, 1);
    ft.savePreviousTimeStep(timeStep, 1);

    ft = surfaceTensionForce_.compute();
    sg = fv::gravity(rho, g_);

    computeRho();
    computeMu();

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) + ib_.eqns(u)
             == fv::laplacian(mu, u) - fv::source(gradP - sg - ft));

    Scalar error = uEqn_.solve();

    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);
    interpolateFaces(fv::INVERSE_VOLUME, gamma);

    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::cn(u, gradGamma, surfaceTensionForce_.n(), gamma, timeStep) + ib_.eqns(gamma) == 0.);

    return gammaEqn_.solve();
}

void FractionalStepMultiphase::computeFaceVelocities(Scalar timeStep)
{
    FractionalStep::computeFaceVelocities(timeStep);

    for(const Face &face: u.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();
        const Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());

        u(face) += timeStep*(ft(face)/rho(face) - g*ft(lCell)/rho(lCell) - (1. - g)*ft(rCell)/rho(rCell));
    }

    for(const Face &face: u.grid.boundaryFaces())
    {
        const Cell& cell = face.lCell();

        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::SYMMETRY:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) += timeStep*(ft(face)/rho(face) - ft(cell)/rho(cell));
            break;
        }
    }
}

void FractionalStepMultiphase::computeMassSource(Scalar timeStep)
{
    FractionalStep::computeMassSource(timeStep);
    return;

    for(const Cell &cell: grid_.fluidCells())
    {
        for(const InteriorLink &nb: cell.neighbours())
            divUStar(cell) -= dot(ft(nb.face()), nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
            divUStar(cell) -= dot(ft(bd.face()), bd.outwardNorm());
    }
}

void FractionalStepMultiphase::computeRho()
{
    const ScalarFiniteVolumeField &w = surfaceTensionForce_.gammaTilde();

    for(const Cell &cell: grid_.activeCells())
        rho(cell) = (1. - w(cell))*rho1_ + w(cell)*rho2_;

    //interpolateFaces(fv::INVERSE_VOLUME, rho);
    harmonicInterpolateFaces(fv::INVERSE_VOLUME, rho);
}

void FractionalStepMultiphase::computeMu()
{
    const ScalarFiniteVolumeField &w = surfaceTensionForce_.gammaTilde();

    for(const Cell &cell: grid_.activeCells())
        mu(cell) = (1. - w(cell))*mu1_ + w(cell)*mu2_;

    interpolateFaces(fv::INVERSE_VOLUME, mu);
    //harmonicInterpolateFaces(mu);
}
