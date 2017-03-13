#include "PisoMultiphase.h"
#include "Cicsam.h"
#include "Plic.h"
#include "Celeste.h"
#include "CrankNicolson.h"
#include "FaceInterpolation.h"
#include "GradientEvaluation.h"
#include "SourceEvaluation.h"
#include "GhostCellImmersedBoundary.h"

PisoMultiphase::PisoMultiphase(const Input &input, const Communicator &comm, FiniteVolumeGrid2D &grid)
        :
        Piso(input, comm, grid),
        gamma(addScalarField(input, "gamma")),
        gradGamma(addVectorField("gradGamma")),
        ft(addVectorField("ft")),
        gradRho(addVectorField("gradRho")),
        sg(addVectorField("sg")),
        gammaEqn_(input, comm, gamma, "gammaEqn")
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

    if (tmp == "CSF")
        surfaceTensionForce_ = std::shared_ptr<SurfaceTensionForce>(new ContinuumSurfaceForce(input, *this));
    else if (tmp == "CELESTE")
        surfaceTensionForce_ = std::shared_ptr<SurfaceTensionForce>(new Celeste(input, *this));
    else
        throw Exception("PisoMultiphase", "PisoMultiphase", "unrecognized surface tension model \"" + tmp + "\".");

    //surfaceTensionForce_->compute();
    computeRho();
    computeMu();

    Scalar sigma = surfaceTensionForce_->sigma();
    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();

    for (const Face &face: grid_.interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_,
                                      sqrt(((rho1_ + rho2_) * delta * delta * delta) / (4 * M_PI * sigma)));
    }

    capillaryTimeStep_ = comm_.min(capillaryTimeStep_);

    comm_.printf("Maximum capillary-wave constrained time-step: %.2e\n", capillaryTimeStep_);
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
            correctPressure();
            correctVelocity();
        }
    }

    solveGammaEqn(timeStep);

    comm_.printf("Max Co = %lf\n", maxCourantNumber(timeStep));

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

    const ScalarFiniteVolumeField &alpha = surfaceTensionForce_->gammaTilde();

    rho.savePreviousTimeStep(0, 1);

    for (const Cell &cell: rho.grid.cells())
    {
        Scalar w = max(0., min(1., alpha(cell)));
        rho(cell) = (1 - w) * rho1_ + w * rho2_;
    }

    grid_.sendMessages(comm_, rho);

    harmonicInterpolateFaces(fv::INVERSE_VOLUME, rho);
    fv::computeInverseWeightedGradient(rho, rho, gradRho);

    for (const Cell &cell: sg.grid.cellZone("fluid"))
        sg(cell) = dot(g_, cell.centroid()) * gradRho(cell);

    for (const Face &face: sg.grid.faces())
        sg(face) = dot(g_, face.centroid()) * gradRho(face);
}

void PisoMultiphase::computeMu()
{
    using namespace std;

    const ScalarFiniteVolumeField &alpha = surfaceTensionForce_->gammaTilde();

    mu.savePreviousTimeStep(0, 1);

    for (const Cell &cell: mu.grid.cells())
    {
        Scalar w = max(0., min(1., alpha(cell)));
        mu(cell) = (1 - w) * mu1_ + w * mu2_;
    }

    grid_.sendMessages(comm_, mu);

    interpolateFaces(fv::INVERSE_VOLUME, mu);
}

Scalar PisoMultiphase::solveUEqn(Scalar timeStep)
{
    ft = surfaceTensionForce_->compute();
    computeRho();
    computeMu();

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho, u, u, 0.5) + ib::gc(ibObjs(), u)
             == cn::laplacian(mu, u, 0.5) + fv::source(ft - gradP - sg));

    Scalar error = uEqn_.solve();

    grid_.sendMessages(comm_, u);

    rhieChowInterpolation();

    return error;
}

Scalar PisoMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);

    switch (interfaceAdvectionMethod_)
    {
        case CICSAM:
            gammaEqn_ = (
                    fv::ddt(gamma, timeStep) + cicsam::cn(u, gradGamma, surfaceTensionForce_->n(), gamma, timeStep) +
                    ib::gc(ibObjs(), gamma) == 0.);
    }

    Scalar error = gammaEqn_.solve();

    grid_.sendMessages(comm_, gamma);

    interpolateFaces(fv::INVERSE_VOLUME, gamma);
    gamma.setBoundaryFaces();

    fv::computeInverseWeightedGradient(rho, gamma, gradGamma);

    grid_.sendMessages(comm_, gradGamma); // Must send gradGamma to other processes for CICSAM to work properly

    return error;
}

void PisoMultiphase::rhieChowInterpolation()
{
    Piso::rhieChowInterpolation();

    for (const Face &face: u.grid.interiorFaces())
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

    for (const Face &face: u.grid.boundaryFaces())
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
