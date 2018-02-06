#include "PisoMultiphase.h"
#include "Cicsam.h"
#include "Celeste.h"
#include "Source.h"

PisoMultiphase::PisoMultiphase(const Input &input,
                               std::shared_ptr<FiniteVolumeGrid2D>& grid)
        :
        Piso(input, grid),
        gamma(addScalarField(input, "gamma")),
        beta(addScalarField("beta")),
        gradGamma(addVectorField(std::make_shared<ScalarGradient>(gamma))),
        gradRho(addVectorField(std::make_shared<ScalarGradient>(rho))),
        sg(addVectorField("sg")),
        gammaEqn_(input, gamma, "gammaEqn")
{
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g", "(0,0)"));

    // volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);

    //- Configuration
    interfaceAdvectionMethod_ = CICSAM;
    const std::string tmp = input.caseInput().get<std::string>("Solver.surfaceTensionModel");

    ft_ = std::make_shared<Celeste>(input, grid_, ib_, rho, mu, u);

    addVectorField(ft_);

    //surfaceTensionForce_->compute();
    computeRho();
    computeMu();

    Scalar sigma = ft_->sigma();
    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();

    for (const Face &face: grid_->interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_,
                                      sqrt(((rho1_ + rho2_) * delta * delta * delta) / (4 * M_PI * sigma)));
    }

    capillaryTimeStep_ = grid_->comm().min(capillaryTimeStep_);

    printf("Maximum capillary-wave constrained time-step: %.2e\n", capillaryTimeStep_);
}

Scalar PisoMultiphase::solve(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    for (size_t innerIter = 0; innerIter < nInnerIterations_; ++innerIter)
    {
        u.savePreviousIteration();
        solveUEqn(timeStep);

        for (size_t pCorrIter = 0; pCorrIter < nPCorrections_; ++pCorrIter)
        {
            solvePCorrEqn();
            correctVelocity();
        }
    }

    solveGammaEqn(timeStep);

    printf("Max Co = %lf\n", maxCourantNumber(timeStep));

    return 0.; // just to get rid of warning
}

Scalar PisoMultiphase::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    return std::min( //- Note: no need to call comm_.min(...) here, args are already global mins
            Piso::computeMaxTimeStep(maxCo, prevTimeStep),
            capillaryTimeStep_
    );
}

//- Protected methods

void PisoMultiphase::computeRho()
{
    using namespace std;

    const ScalarFiniteVolumeField &alpha = gamma;

    rho.savePreviousTimeStep(0, 1);

    for (const Cell &cell: rho.grid()->cells())
    {
        Scalar w = max(0., min(1., alpha(cell)));
        rho(cell) = (1 - w) * rho1_ + w * rho2_;
    }

    grid_->sendMessages(rho);

    rho.interpolateFaces();
    gradRho.compute(fluid_);

    for (const Cell &cell: sg.grid()->cellZone("fluid"))
        sg(cell) = dot(g_, cell.centroid()) * gradRho(cell);

    for (const Face &face: sg.grid()->faces())
        sg(face) = dot(g_, face.centroid()) * gradRho(face);
}

void PisoMultiphase::computeMu()
{
    using namespace std;

    const ScalarFiniteVolumeField &alpha = gamma;

    mu.savePreviousTimeStep(0, 1);

    for (const Cell &cell: mu.grid()->cells())
    {
        Scalar w = max(0., min(1., alpha(cell)));
        mu(cell) = (1 - w) * mu1_ + w * mu2_;
    }

    grid_->sendMessages(mu);
    mu.interpolateFaces();
}

Scalar PisoMultiphase::solveUEqn(Scalar timeStep)
{
    ft_->compute(gamma, gradGamma);
    computeRho();
    computeMu();

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) + ib_->bcs(u)
             == fv::laplacian(mu, u) + src::src(*ft_ - gradP - sg, fluid_));

    Scalar error = uEqn_.solve();

    grid_->sendMessages(u);

    rhieChowInterpolation();

    return error;
}

Scalar PisoMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);

    cicsam::beta(u, gradGamma, gamma, timeStep, beta, 0.5);

    switch (interfaceAdvectionMethod_)
    {
        case CICSAM:
            gammaEqn_ = (
                    fv::ddt(gamma, timeStep) + cicsam::div(u, beta, gamma, 0.5) +
                    ib_->bcs(gamma) == 0.);
    }

    Scalar error = gammaEqn_.solve();

    grid_->sendMessages(gamma);

    gamma.interpolateFaces([this](const Face& face){
        return dot(u(face), face.outwardNorm(face.lCell().centroid())) > 0 ? 1. - beta(face): beta(face);
    });

    gradGamma.compute(fluid_);

    grid_->sendMessages(gradGamma); // Must send gradGamma to other processes for CICSAM to work properly

    return error;
}

void PisoMultiphase::rhieChowInterpolation()
{
    Piso::rhieChowInterpolation();
    const auto& ft = *ft_;

    for (const Face &face: u.grid()->interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();

        const Scalar df = d(face);
        const Scalar g = rCell.volume() / (lCell.volume() + rCell.volume());

        //        if(!(grid_.cellZone("fluid").isInGroup(lCell) && grid_.cellZone("fluid").isInGroup(rCell)))
        //            continue;

        u(face) += df * ft(face) - (g * d(lCell) * ft(lCell) + (1. - g) * d(rCell) * ft(rCell))
                   + df * sg(face) - (g * d(lCell) * sg(lCell) + (1. - g) * d(rCell) * sg(rCell));
    }

    for (const Face &face: u.grid()->boundaryFaces())
    {
        const Cell &cellP = face.lCell();
        const Scalar df = d(face);

        switch (u.boundaryType(face))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                u(face) += df * ft(face) - d(cellP) * ft(cellP)
                           + df * sg(face) - d(cellP) * sg(cellP);
                break;
        }
    }
}

//- External functions
