#include "FractionalStepMultiphase.h"
#include "Cicsam.h"
#include "CrankNicolson.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"
#include "SourceEvaluation.h"

FractionalStepMultiphase::FractionalStepMultiphase(const Input &input, const Communicator& comm, FiniteVolumeGrid2D& grid)
    :
      FractionalStep(input, comm, grid),
      gamma(addScalarField(input, "gamma")),
      gradGamma(addVectorField("gradGamma")),
      ft(addVectorField("ft")),
      sg(addVectorField("sg")),
      gradRho(addVectorField("gradRho")),
      gammaEqn_(input, comm, gamma, "gammaEqn"),
      surfaceTensionForce_(input, *this)
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g"));

    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);

    surfaceTensionForce_.compute();
    computeRho();
    computeMu();

    Scalar sigma = surfaceTensionForce_.sigma();
    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();

    for(const Face& face: grid_.interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_, sqrt(((rho1_ + rho2_)*delta*delta*delta)/(4*M_PI*sigma)));
    }

    printf("Maximum capillary-wave constrained time-step: %.2e\n", capillaryTimeStep_);
}

Scalar FractionalStepMultiphase::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    solveGammaEqn(timeStep);

    printf("Max Co = %lf\n", maxCourantNumber(timeStep));

    return 0.;
}

Scalar FractionalStepMultiphase::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    return std::min(
                FractionalStep::computeMaxTimeStep(maxCo, prevTimeStep),
                capillaryTimeStep_
                );
}

//- Protected methods

Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    ft.savePreviousTimeStep(timeStep, 1);
    ft = surfaceTensionForce_.compute();

    computeRho();
    computeMu();

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho.prev(0)*u, u) + ib_.eqns(u)
             == fv::laplacian(mu, u) - fv::source(gradP - ft - sg.prev(0)));

    Scalar error = uEqn_.solve();

    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);

    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::cn(u, gradGamma, surfaceTensionForce_.n(), gamma, timeStep) + ib_.eqns(gamma) == 0.);

    gamma.setBoundaryFaces();
    fv::computeInverseWeightedGradient(rho, gamma, gradGamma);

    return gammaEqn_.solve();
}

void FractionalStepMultiphase::computeFaceVelocities(Scalar timeStep)
{
    FractionalStep::computeFaceVelocities(timeStep);

    const ScalarFiniteVolumeField& rho0 = rho.prev(0);
    const VectorFiniteVolumeField& ft0 = ft.prev(0);
    const VectorFiniteVolumeField& sg0 = sg.prev(0);

    for(const Face &face: u.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        if(!(lCell.isFluidCell() && rCell.isFluidCell()))
            continue;

        const Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());

        u(face) += timeStep*(ft(face)/rho(face) - g*ft0(lCell)/rho0(lCell) - (1. - g)*ft0(rCell)/rho0(rCell)
                             + sg(face)/rho(face) - g*sg0(lCell)/rho0(lCell) - (1. - g)*sg0(rCell)/rho0(rCell)
                             );
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
            u(face) += timeStep*(ft(face)/rho(face) - ft(cell)/rho0(cell)
                                 + sg(face)/rho(face) - ft(cell)/rho0(cell)
                                 );
            break;
        }
    }
}

void FractionalStepMultiphase::correctVelocity(Scalar timeStep)
{
    const ScalarFiniteVolumeField& rho0 = rho.prev(0);
    const VectorFiniteVolumeField& gradP0 = gradP.prev(0);
    const VectorFiniteVolumeField& ft0 = ft.prev(0);
    const VectorFiniteVolumeField& sg0 = sg.prev(0);

    for(const Cell &cell: grid_.fluidCells())
        u(cell) -= timeStep*(gradP(cell)/rho(cell) - gradP0(cell)/rho0(cell)
                             - ft(cell)/rho(cell) + ft0(cell)/rho0(cell)
                             - sg(cell)/rho(cell) + sg0(cell)/rho0(cell)
                             );

    for(const Face &face: grid_.interiorFaces())
        u(face) -= timeStep/rho(face)*gradDp(face);

    for(const Face &face: grid_.boundaryFaces())
    {
        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall)*nWall/nWall.magSqr();
        }
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) -= timeStep/rho(face)*gradDp(face);
            break;
        };
    }
}

void FractionalStepMultiphase::computeRho()
{
    const ScalarFiniteVolumeField &w = surfaceTensionForce_.gammaTilde();

    rho.savePreviousTimeStep(0, 1);

    for(const Cell &cell: grid_.activeCells())
    {
        Scalar g = std::max(0., std::min(1., w(cell)));
        rho(cell) = (1. - g)*rho1_ + g*rho2_;
    }

    harmonicInterpolateFaces(fv::INVERSE_VOLUME, rho);

    fv::computeInverseWeightedGradient(rho, rho, gradRho);

    sg.savePreviousTimeStep(0, 1);

    for(const Cell& cell: sg.grid.fluidCells())
        sg(cell) = dot(g_, -cell.centroid())*gradRho(cell);

    for(const Face& face: sg.grid.faces())
        sg(face) = dot(g_, -face.centroid())*gradRho(face);
}

void FractionalStepMultiphase::computeMu()
{
    const ScalarFiniteVolumeField &w = surfaceTensionForce_.gammaTilde();

    mu.savePreviousTimeStep(0, 1);

    for(const Cell &cell: grid_.activeCells())
    {
        Scalar g = std::max(0., std::min(1., w(cell)));
        mu(cell) = (1. - g)*mu1_ + g*mu2_;
    }

    interpolateFaces(fv::INVERSE_VOLUME, mu);
}
