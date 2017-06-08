#include "FractionalStepMultiphase.h"
#include "CrankNicolson.h"
#include "Cicsam.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"
#include "SourceEvaluation.h"
#include "Source.h"

FractionalStepMultiphase::FractionalStepMultiphase(const Input &input, const Communicator &comm,
                                                   FiniteVolumeGrid2D &grid)
        :
        FractionalStep(input, comm, grid),
        gamma(addScalarField(input, "gamma")),
        gammaFlux(addScalarField("gammaFlux")),
        gradGamma(addVectorField("gradGamma")),
        ft(addVectorField("ft")),
        sg(addVectorField("sg")),
        gradRho(addVectorField("gradRho")),
        rhoU(addVectorField("rhoU")),
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

    for (const Face &face: grid_.interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_,
                                      sqrt(((rho1_ + rho2_) * delta * delta * delta) / (4 * M_PI * sigma)));
    }

    capillaryTimeStep_ = comm.min(capillaryTimeStep_);
    comm.printf("CICSAM blending constant (k): %.2f\n", cicsamBlending_);
    comm.printf("Maximum capillary-wave constrained time-step: %.2e\n", capillaryTimeStep_);
}

void FractionalStepMultiphase::initialize()
{
    //- Ensure the computation starts with a valid gamma field
    fv::computeGradient(fv::FACE_TO_CELL, fluid_, gamma, gradGamma);
    ft = surfaceTensionForce_.compute();
    ft.savePreviousTimeStep(0, 1);

    //- Be careful about this, mu and rho must have valid time histories!!!
    computeRho();
    computeMu();
    computeRho();
    computeMu();
}

Scalar FractionalStepMultiphase::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep); //- Solve a sharp gamma equation
    solveUEqn(timeStep); //- Solve a momentum prediction

    solvePEqn(timeStep); //- Solve the pressure equation, using sharp value of rho
    correctVelocity(timeStep);

    comm_.printf("Max Co = %lf\n", maxCourantNumber(timeStep));
    comm_.printf("Max absolute velocity divergence error = %.4e\n", maxDivergenceError());

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
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u) + ib_.bcs(u)
             == cn::laplacian(mu, u, 0.5) - fv::source(gradP - sg.oldField(0) - ft.oldField(0)));

    Scalar error = uEqn_.solve();
    grid_.sendMessages(comm_, u);
    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStepMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho, p) + ib_.bcs(p) ==
             source::div(u));

    Scalar error = pEqn_.solve();
    grid_.sendMessages(comm_, p);

    //- Compute pressure gradient
    gradP.savePreviousTimeStep(timeStep, 1);

    //- Weighted gradients greatly reduce the effect of large pressure differences
    fv::computeInverseWeightedGradient(rho, p, gradP);
    grid_.sendMessages(comm_, gradP);

    return error;
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);

    cicsam::interpolateFaces(u, gradGamma, gamma, timeStep, cicsamBlending_);
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, gamma) +
                 ib_.bcs(gamma) == 0.);

    Scalar error = gammaEqn_.solve();

    //- While this may affect mass conservation, it prevents issues at high density ratios
    for (const Cell &cell: grid_.localActiveCells())
        gamma(cell) = std::max(std::min(gamma(cell), 1.), 0.);

    grid_.sendMessages(comm_, gamma);

    for (const Face &face: grid_.faces())
        rhoU(face) = ((1. - gamma(face)) * rho1_ + gamma(face) * rho2_) * u(face);

    //- Recompute gradGamma for next cicsam iteration
    fv::computeGradient(fv::FACE_TO_CELL, fluid_, gamma, gradGamma);
    grid_.sendMessages(comm_, gradGamma); //- In case donor cell is on another proc

    // compute surface tension for at new time level
    ft.savePreviousTimeStep(timeStep, 1);
    ft = surfaceTensionForce_.compute();

    computeRho();
    computeMu();

    return error;
}

void FractionalStepMultiphase::computeFaceVelocities(Scalar timeStep)
{
    const ScalarFiniteVolumeField &rho0 = rho.oldField(0);
    const VectorFiniteVolumeField &sg0 = sg.oldField(0);
    const VectorFiniteVolumeField &ft0 = ft.oldField(0);

    for (const Face &face: u.grid.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();

        Scalar l1 = dot(face.centroid() - lCell.centroid(), (rCell.centroid() - lCell.centroid()).unitVec());
        Scalar l2 = dot(rCell.centroid() - face.centroid(), (rCell.centroid() - lCell.centroid()).unitVec());
        Scalar g = mu(lCell) / l1 / (mu(lCell) / l1 + mu(rCell) / l2);

        u(face) = g * (u(lCell) + timeStep / rho0(lCell) * (gradP(lCell) - sg0(lCell) - ft0(lCell)))
                  + (1. - g) * (u(rCell) + timeStep / rho0(rCell) * (gradP(rCell) - sg0(rCell) - ft0(rCell)))
                  + timeStep / rho(face) * (sg(face) + ft(face));
    }

    for (const Face &face: u.grid.boundaryFaces())
    {
        const Cell &cell = face.lCell();

        switch (u.boundaryType(face))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                u(face) = u(cell) + timeStep / rho0(cell) * (gradP(cell) - sg0(cell) - ft0(cell)) +
                          timeStep / rho(face) * (sg(face) + ft(face));
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = face.outwardNorm(cell.centroid());
                u(face) = u(cell) - dot(u(cell), nWall) * nWall / nWall.magSqr();
            }
                break;

                //default:
                //throw Exception("FractionalStep", "computeFaceVelocities", "unrecognized boundary condition type.");
        };
    }
}

void FractionalStepMultiphase::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField &gradP0 = gradP.oldField(0);
    const VectorFiniteVolumeField &sg0 = sg.oldField(0);
    const VectorFiniteVolumeField &ft0 = ft.oldField(0);
    const ScalarFiniteVolumeField &rho0 = rho.oldField(0);

    for (const Cell &cell: fluid_)
        u(cell) -= timeStep * ((gradP(cell) - sg(cell) - ft(cell)) / rho(cell) -
                               (gradP0(cell) - sg0(cell) - ft0(cell)) / rho0(cell));

    grid_.sendMessages(comm_, u);

    for (const Face &face: grid_.interiorFaces())
        u(face) -= timeStep / rho(face) * gradP(face);

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
                u(face) -= timeStep / rho(face) * gradP(face);
                break;
        };
    }
}

void FractionalStepMultiphase::computeRho()
{
    rho.savePreviousTimeStep(0, 1);

    //- Update density
    for (const Cell &cell: grid_.cells())
    {
        Scalar g = std::max(0., std::min(1., gamma(cell)));
        rho(cell) = (1. - g) * rho1_ + g * rho2_;
    }

    grid_.sendMessages(comm_, rho);

    for (const Face &face: grid_.boundaryFaces())
    {
        Scalar g = std::max(0., std::min(1., gamma(face)));
        rho(face) = (1. - g) * rho1_ + g * rho2_;
    }

    //- Recompute rho_f to ensure consistency with momentum advection
    fv::interpolateFaces(fv::INVERSE_VOLUME, rho);
    fv::computeInverseWeightedGradient(rho, rho, gradRho);

    //- Update the gravitational source term
    sg.savePreviousTimeStep(0., 1);
    for (const Cell &cell: fluid_)
        sg(cell) = dot(g_, -cell.centroid()) * gradRho(cell);

    grid_.sendMessages(comm_, sg);

    for (const Face &face: sg.grid.faces())
        sg(face) = dot(g_, -face.centroid()) * gradRho(face);
}

void FractionalStepMultiphase::computeMu()
{
    mu.savePreviousTimeStep(0, 1);

    //- Update viscosity
    for (const Cell &cell: grid_.cells())
    {
        Scalar g = std::max(0., std::min(1., gamma(cell)));
        mu(cell) = (1. - g) * mu1_ + g * mu2_;
    }

    grid_.sendMessages(comm_, mu);
    fv::interpolateFaces(fv::INVERSE_VOLUME, mu);
}

void FractionalStepMultiphase::checkMassFluxConsistency(Scalar timeStep)
{
    Scalar max = 0.;

    for (const Face &face: grid_.faces())
    {
        Scalar gammaMom = (rho(face) - rho1_) / (rho2_ - rho1_);
        Scalar gammaVof = gamma(face);

        Scalar error = fabs((gammaMom - gammaVof));

        if (error > max)
            max = error;
    }

    comm_.printf("Max mass flux error = %.2lf.\n", max);
}

void FractionalStepMultiphase::computeGammaFux(Scalar timeStep)
{
    gammaFlux.fill(0.);

    for (const Cell &cell: fluid_)
    {
        for (const InteriorLink &nb: cell.neighbours())
        {
            gammaFlux(cell) -= gamma(nb.face()) * dot(u(nb.face()), nb.outwardNorm()) * timeStep / cell.volume();
        }
    }
}
