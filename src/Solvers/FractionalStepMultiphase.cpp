#include "FractionalStepMultiphase.h"
#include "Cicsam.h"
#include "CrankNicolson.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"
#include "SourceEvaluation.h"
#include "GhostCellImmersedBoundary.h"
#include "TimeDerivative.h"

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
    cicsamBlending_ = input.caseInput().get<Scalar>("Solver.cicsamBlending", 1.);
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1");
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2");
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1");
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2");

    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g"));

    //volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);
    Scalar sigma = surfaceTensionForce_.sigma();
    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();

    for(const Face& face: grid_.interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_, sqrt(((rho1_ + rho2_)*delta*delta*delta)/(4*M_PI*sigma)));
    }

    capillaryTimeStep_ = comm.min(capillaryTimeStep_);
    comm.printf("CICSAM blending constant (k): %.2f\n", cicsamBlending_);
    comm.printf("Maximum capillary-wave constrained time-step: %.2e\n", capillaryTimeStep_);
}

void FractionalStepMultiphase::initialize()
{
    surfaceTensionForce_.compute();
    computeRho();
    computeMu();
    computeRho();
    computeMu();
    u.savePreviousTimeStep(0, 1);
}

Scalar FractionalStepMultiphase::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep);

    ft.savePreviousTimeStep(timeStep, 1);
    ft = surfaceTensionForce_.compute();

    //- Update the gravitational source term
    sg.savePreviousTimeStep(timeStep, 1);
    for(const Cell& cell: sg.grid.cellZone("fluid"))
        sg(cell) = dot(g_, -cell.centroid())*gradRho(cell);

    for(const Face& face: sg.grid.faces())
        sg(face) = dot(g_, -face.centroid())*gradRho(face);

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    comm_.printf("Max Co = %lf\n", maxCourantNumber(timeStep));

    return 0.;
}

Scalar FractionalStepMultiphase::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    return std::min( //- Both args have already been globalling minimized
                     FractionalStep::computeMaxTimeStep(maxCo, prevTimeStep),
                     capillaryTimeStep_
                     );
}

//- Protected methods
Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho, u, u, 0.) + ib::gc(ibObjs(), u)
             == cn::laplacian(mu, u, 1.5) - fv::source(gradP - ft.prev(0) - sg.prev(0)));

    Scalar error = uEqn_.solve();
    grid_.sendMessages(comm_, u);
    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::cn(u, gradGamma, surfaceTensionForce_.n(), gamma, timeStep, 0.5, cicsamBlending_) + ib::gc(ibObjs(), gamma) == 0.);
    //gammaEqn_ = (fv::ddt(gamma, timeStep) + fv::div(u, gamma) + ib::gc(ibObjs(), gamma) == 0.);

    Scalar error = gammaEqn_.solve();
    grid_.sendMessages(comm_, gamma);

    interpolateFaces(fv::INVERSE_VOLUME, gamma);
    //gamma.setBoundaryFaces();

    //cicsam::interpolateFaces(u, gradGamma, surfaceTensionForce_.n(), gamma, timeStep, cicsamBlending_);
    //gamma.setBoundaryFaces();

    computeRho();
    computeMu();

    fv::computeInverseWeightedGradient(rho, gamma, gradGamma);
    grid_.sendMessages(comm_, gradGamma); // Must send gradGamma to other processes for CICSAM to work properly

    return error;
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
        const Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());

        u(face) += timeStep*(ft0(face)/rho(face) - g*ft0(lCell)/rho(lCell) - (1. - g)*ft0(rCell)/rho(rCell)
                             + sg0(face)/rho(face) - g*sg0(lCell)/rho(lCell) - (1. - g)*sg0(rCell)/rho(rCell)
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
            u(face) += timeStep*(ft0(face)/rho(face) - ft0(cell)/rho(cell)
                                 + sg0(face)/rho(face) - sg0(cell)/rho(cell)
                                 );
            break;
        }
    }
}

void FractionalStepMultiphase::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField& gradP0 = gradP.prev(0);
    const VectorFiniteVolumeField& ft0 = ft.prev(0);
    const VectorFiniteVolumeField& sg0 = sg.prev(0);

    for (const Cell &cell: grid_.cellZone("fluid"))
        u(cell) -= timeStep / rho(cell) * (gradP(cell) - gradP0(cell) - ft(cell) + ft0(cell) - sg(cell) + sg0(cell));

    grid_.sendMessages(comm_, u);

    for (const Face &face: grid_.interiorFaces())
        u(face) -= timeStep / rho(face) * (gradP(face) - gradP0(face) - ft(face) + ft0(face) - sg(face) + sg0(face));

    for (const Face &face: grid_.boundaryFaces())
    {
        switch (u.boundaryType(face))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
                u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
            }
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                u(face) -= timeStep / rho(face) * (gradP(face) - gradP0(face) - ft(face) + ft0(face) - sg(face) + sg0(face));
                break;
        };
    }
}

void FractionalStepMultiphase::computeRho()
{
    const ScalarFiniteVolumeField &w = gamma;
    rho.savePreviousTimeStep(0, 1);

    //- Update density
    for(const Cell &cell: grid_.localActiveCells())
    {
        Scalar g = std::max(0., std::min(1., w(cell)));
        rho(cell) = (1. - g)*rho1_ + g*rho2_;
    }

    grid_.sendMessages(comm_, rho);
    //harmonicInterpolateFaces(fv::INVERSE_VOLUME, rho);
    //interpolateFaces(fv::INVERSE_VOLUME, rho);

    for(const Face& face: grid_.faces())
    {
        Scalar g = std::max(0., std::min(1., w(face)));
        rho(face) = (1. - g)*rho1_ + g*rho2_;
    }

    fv::computeInverseWeightedGradient(rho, rho, gradRho);
}

void FractionalStepMultiphase::computeMu()
{
    const ScalarFiniteVolumeField &w = gamma;

    mu.savePreviousTimeStep(0, 1);

    //- Update viscosity
    for(const Cell &cell: grid_.localActiveCells())
    {
        Scalar g = std::max(0., std::min(1., w(cell)));
        mu(cell) = (1. - g)*mu1_ + g*mu2_;
    }

    grid_.sendMessages(comm_, mu);

    for(const Face& face: grid_.faces())
    {
        Scalar g = std::max(0., std::min(1., w(face)));
        mu(face) = (1. - g)*mu1_ + g*mu2_;
    }
}
