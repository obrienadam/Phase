#include "FractionalStepMultiphase.h"
#include "Cicsam.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

FractionalStepMultiphase::FractionalStepMultiphase(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      FractionalStep(grid, input),
      gamma(addScalarField(input, "gamma")),
      gradGamma(addVectorField("gradGamma")),
      ft(addVectorField("ft")),
      gammaEqn_(gamma, "gammaEqn"),
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

    for(const ForceIntegrator &fi: forceIntegrators_)
        fi.integrate();

    for(const VolumeIntegrator &vi: volumeIntegrators_)
        vi.integrate();

    printf("Max Co = %lf\n", courantNumber(timeStep));

    return 0.;
}

//- Protected methods

Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    ft = surfaceTensionForce_.compute();
    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(u, timeStep) + cn::div(u, u) + ib_.eqns(u)
             == ab::laplacian(mu/rho, u) - fv::source(gradP/rho) + fv::source(sg/rho) + fv::source(ft/rho));

    Scalar error = uEqn_.solve();

    for(const Cell& cell: u.grid.fluidCells())
        u(cell) += timeStep*(gradP(cell) - sg(cell))/rho(cell);

    balancedForceInterpolation(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);
    interpolateFaces(fv::INVERSE_VOLUME, gamma);

    gammaEqn_ = (cicsam::div(u, gradGamma, ib_.ibObjs(), gamma, timeStep, cicsam::HC) == 0.);

    Scalar error = gammaEqn_.solve();

    computeRho();
    computeMu();

    return error;
}

void FractionalStepMultiphase::balancedForceInterpolation(Scalar timeStep)
{
    const Scalar sigma = surfaceTensionForce_.sigma();

    for(const Face &face: u.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        const Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());
        const Scalar kf = surfaceTensionForce_.kappa()(face);
        const Vector2D rc = rCell.centroid() - lCell.centroid();

        u(face) = g*(u(lCell) - timeStep*ft(lCell)/rho(lCell))
                + (1. - g)*(u(rCell) - timeStep*ft(rCell)/rho(rCell))
                + timeStep/rho(face)*sigma*kf*(gamma(rCell) - gamma(lCell))*rc/dot(rc, rc);
    }

    for(const Face &face: u.grid.boundaryFaces())
    {
        const Cell& cell = face.lCell();

        const Scalar kf = surfaceTensionForce_.kappa()(face);
        const Vector2D rf = face.centroid() - cell.centroid();

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            u(face) = u(cell) - dot(u(cell), nWall)*nWall/nWall.magSqr();
            break;
        }

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) = u(cell) - timeStep*ft(cell)/rho(cell)
                    + timeStep/rho(face)*sigma*kf*(gamma(face) - gamma(cell))*rf/dot(rf, rf);
            break;
        }
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
