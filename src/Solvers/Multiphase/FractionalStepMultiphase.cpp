#include "FractionalStepMultiphase.h"
#include "Cicsam.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"
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

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) + ib_.eqns(u)
             == fv::laplacian(mu, u) - fv::source(gradP) + fv::source(sg) + fv::source(ft));

    Scalar error = uEqn_.solve();
    computeAdvectingVelocity(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);

    interpolateFaces(fv::INVERSE_VOLUME, gamma);

    //gammaEqn_ = (cicsam::div(u, gamma, timeStep, cicsam::HC) + ib_.eqns(gamma) == 0.);
    gammaEqn_ = (cicsam::div(u, gradGamma, ib_.ibObjs(), gamma, timeStep, cicsam::HC) == 0.);

    Scalar error = gammaEqn_.solve();

    computeRho();
    computeMu();

    return error;
}

void FractionalStepMultiphase::computeAdvectingVelocity(Scalar timeStep)
{
    //FractionalStep::computeAdvectingVelocity(timeStep);
    const Scalar sigma = surfaceTensionForce_.sigma();

    for(const Face &face: u.grid.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();

        Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());
        Vector2D rc = rCell.centroid() - lCell.centroid();

        const Scalar rhof = rho(face);
        const Scalar kf = surfaceTensionForce_.kappa()(face);
        const Scalar p0 = p(lCell), p1 = p(rCell);
        const Scalar g0 = gamma(lCell), g1 = gamma(rCell);

        u(face) = g*(u(lCell) + timeStep*(gradP(lCell) - sg(lCell) - ft(lCell))/rho(lCell))
                + (1. - g)*(u(rCell) + timeStep*(gradP(rCell) - sg(rCell) - ft(rCell))/rho(rCell))
                - timeStep/rhof*(p1 - p0)*rc/rc.magSqr()
                + timeStep/rhof*g_;
                + timeStep/rhof*sigma*kf*(g1 - g0)*rc/rc.magSqr();
    }

    for(const Face &face: u.grid.boundaryFaces())
    {
        const Cell &cell = face.lCell();

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT: case VectorFiniteVolumeField::OUTFLOW:
        {
            Vector2D rf = face.centroid() - cell.centroid();
            const Scalar rhof = rho(face);
            const Scalar kf = surfaceTensionForce_.kappa()(face);
            const Scalar p0 = p(cell), p1 = p(face);
            const Scalar g0 = gamma(cell), g1 = gamma(face);

            u(face) = u(cell)
                    + timeStep*(gradP(cell) - sg(cell) - ft(cell))/rho(cell)
                    - timeStep/rhof*(p1 - p0)*rf/rf.magSqr()
                    + timeStep/rhof*g_
                    + timeStep/rhof*sigma*kf*(g1 - g0)*rf/rf.magSqr();
        }
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            u(face) = u[face.lCell().id()] - dot(u[face.lCell().id()], nWall)*nWall/nWall.magSqr();
            break;
        }

        default:
            throw Exception("FractionalStepMultiphase", "computeAdvectingVelocity", "unrecongnized boundary condition type.");
        }
    }
}

void FractionalStepMultiphase::computeRho()
{
    const ScalarFiniteVolumeField &w = surfaceTensionForce_.gammaTilde();

    for(const Cell &cell: grid_.activeCells())
        rho(cell) = (1. - w(cell))*rho1_ + w(cell)*rho2_;

    interpolateFaces(fv::INVERSE_VOLUME, rho);
    //harmonicInterpolateFaces(rho);
}

void FractionalStepMultiphase::computeMu()
{
    const ScalarFiniteVolumeField &w = surfaceTensionForce_.gammaTilde();

    for(const Cell &cell: grid_.activeCells())
        mu(cell) = (1. - w(cell))*mu1_ + w(cell)*mu2_;

    interpolateFaces(fv::INVERSE_VOLUME, mu);
    //harmonicInterpolateFaces(mu);
}
