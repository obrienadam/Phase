#include "Algorithm.h"
#include "FractionalStepIncrementalMultiphase.h"
#include "Cicsam.h"
#include "Hric.h"
#include "FaceInterpolation.h"
#include "Source.h"
#include "GhostCellImmersedBoundaryObjectContactLineTracker.h"

FractionalStepIncrementalMultiphase::FractionalStepIncrementalMultiphase(const Input &input,
                                                                         std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStepIncremental(input, grid),
        gamma(addScalarField(input, "gamma")),
        rho(addScalarField("rho")),
        mu(addScalarField("mu")),
        gradGamma(addVectorField(std::make_shared<ScalarGradient>(gamma))),
        gradRho(addVectorField(std::make_shared<ScalarGradient>(rho))),
        ft(addVectorField(std::make_shared<Celeste>(input, ib_, gamma, rho, mu, u, gradGamma))),
        sg(addVectorField("sg")),
        rhoU(addVectorField("rhoU")),
        gammaEqn_(input, gamma, "gammaEqn")
{
    cicsamBlending_ = input.caseInput().get<Scalar>("Solver.cicsamBlending", 1.);
    rho1_ = input.caseInput().get<Scalar>("Properties.rho1", rho_);
    rho2_ = input.caseInput().get<Scalar>("Properties.rho2", rho_);
    mu1_ = input.caseInput().get<Scalar>("Properties.mu1", mu_);
    mu2_ = input.caseInput().get<Scalar>("Properties.mu2", mu_);
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g", "(0,0)"));

    addScalarField(ft.kappaPtr());
    addScalarField(ft.gammaTildePtr());
    addVectorField(ft.nPtr());
    addVectorField(ft.gradGammaTildePtr());

    capillaryTimeStep_ = std::numeric_limits<Scalar>::infinity();
    for (const Face &face: grid_->interiorFaces())
    {
        Scalar delta = (face.rCell().centroid() - face.lCell().centroid()).mag();
        capillaryTimeStep_ = std::min(capillaryTimeStep_,
                                      sqrt(((rho1_ + rho2_) * pow(delta, 3)) / (4 * M_PI * ft.sigma())));
    }

    capillaryTimeStep_ = grid_->comm().min(capillaryTimeStep_);
    printf("CICSAM blending constant (k): %.2f\n", cicsamBlending_);
    printf("Maximum capillary-wave constrained time-step: %.2e\n", capillaryTimeStep_);
}

void FractionalStepIncrementalMultiphase::initialize()
{
    FractionalStepIncremental::initialize();

    //- Ensure the computation starts with a valid gamma field
    gradGamma.compute(fluid_);

    //- Make sure a time history exists
    updateProperties(0.);
    updateProperties(0.);
}

Scalar FractionalStepIncrementalMultiphase::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep); //- Solve a sharp gamma equation
    solveUEqn(timeStep); //- Solve a momentum prediction
    solvePEqn(timeStep); //- Solve the pressure equation, using sharp value of rho
    correctVelocity(timeStep);

    grid_->comm().printf("Max Co = %lf\n", maxCourantNumber(timeStep));
    grid_->comm().printf("Max absolute velocity divergence error = %.4e\n", maxDivergenceError());

    return 0.;
}

Scalar FractionalStepIncrementalMultiphase::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    return std::min( //- Both args have already been globally minimized
            FractionalStepIncremental::computeMaxTimeStep(maxCo, prevTimeStep),
            capillaryTimeStep_
    );
}

Scalar FractionalStepIncrementalMultiphase::solveGammaEqn(Scalar timeStep)
{
    auto beta = cicsam::beta(u, gradGamma, gamma, timeStep);

    gamma.savePreviousTimeStep(timeStep, 1);
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, beta, gamma, 0.5) + ib_.contactLineBcs(ft, gamma) == 0.);

    //- Solve and update faces
    Scalar error = gammaEqn_.solve();
    grid_->sendMessages(gamma);
    gamma.interpolateFaces();

    //- Update gradient
    gradGamma.compute(fluid_);
    grid_->sendMessages(gradGamma); //- In case donor cell is on another proc

    //- Update mass flux
    rhoU.savePreviousTimeStep(timeStep, 1);
    rhoU.oldField(0).computeInteriorFaces([this, &beta](const Face &f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar g = (1. - beta(f)) * gamma.oldField(0)(d) + beta(f) * gamma.oldField(0)(a);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    rhoU.oldField(0).computeBoundaryFaces([this](const Face &f) {
        Scalar g = gamma.oldField(0)(f);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    rhoU.computeInteriorFaces([this, &beta](const Face &f) {
        Scalar flux = dot(u(f), f.outwardNorm());
        const Cell &d = flux > 0. ? f.lCell() : f.rCell();
        const Cell &a = flux > 0. ? f.rCell() : f.lCell();
        Scalar g = (1. - beta(f)) * gamma(d) + beta(f) * gamma(a);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    rhoU.computeBoundaryFaces([this](const Face &f) {
        Scalar g = gamma(f);
        return ((1. - g) * rho1_ + g * rho2_) * u(f);
    });

    updateProperties(timeStep);

    return error;
}

//- Protected methods
Scalar FractionalStepIncrementalMultiphase::solveUEqn(Scalar timeStep)
{
    const VectorFiniteVolumeField &sg0 = sg.oldField(0);
    const VectorFiniteVolumeField &ft0 = ft.oldField(0);

    gradP.faceToCell(rho, rho.oldField(0), fluid_);

    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u, 0.5) + ib_.solidVelocity(u)
             == fv::laplacian(mu, u, 0.5) - src::src(gradP - sg0 - ft0, fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    for (const Face &f: grid_->interiorFaces())
    {
        Scalar g = f.volumeWeight();
        const Cell &l = f.lCell();
        const Cell &r = f.rCell();

        u(f) = g * (u(l) + timeStep / rho(l) * (gradP(l) - sg0(l) - ft0(l)))
               + (1. - g) * (u(r) + timeStep / rho(r) * (gradP(r) - sg0(r) - ft0(r)))
               - timeStep / rho(f) * (gradP(f) - sg(f) - ft(f));
    }

    for (const Patch &patch: u.grid().patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &f: patch)
                {
                    const Cell &l = f.lCell();
                    u(f) = u(l) + timeStep / rho(l) * (gradP(l) - sg0(l) - ft0(l))
                           - timeStep / rho(f) * (gradP(f) - sg(f) - ft(f));
                }
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &f: patch)
                    u(f) = u(f.lCell()) - dot(u(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
                break;
        }

    return error;
}

Scalar FractionalStepIncrementalMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho, p) + ib_.bcs(p) == src::div(u) + src::laplacian(timeStep / rho, p));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    p.setBoundaryFaces();

    gradP.savePreviousTimeStep(timeStep, 1);
    gradP.computeFaces();
    gradP.faceToCell(rho, rho, fluid_);

    return error;
}

void FractionalStepIncrementalMultiphase::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField &gradP0 = gradP.oldField(0);
    VectorFiniteVolumeField &sg0 = sg.oldField(0);
    VectorFiniteVolumeField &ft0 = ft.oldField(0);

    for (const Cell &cell: fluid_)
        u(cell) -= timeStep / rho(cell) * (gradP(cell) - gradP0(cell)
                                           - sg(cell) + sg0(cell)
                                           - ft(cell) + ft0(cell));

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho(face) * (gradP(face) - gradP0(face));

    for (const Patch &patch: grid_->patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &face: patch)
                    u(face) -= timeStep / rho(face) * (gradP(face) - gradP0(face));
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &f: patch)
                    u(f) = u(f.lCell()) - dot(u(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
                break;
        }

    grid_->sendMessages(u);
}

void FractionalStepIncrementalMultiphase::updateProperties(Scalar timeStep)
{
    //- Update density
    rho.savePreviousTimeStep(timeStep, 1);
    rho.computeCells([this](const Cell &cell) {
        Scalar g = gamma(cell);
        return (1. - g) * rho1_ + g * rho2_;
    });

    rho.computeFaces([this](const Face &face) {
        Scalar g = gamma(face);
        return (1. - g) * rho1_ + g * rho2_;
    });

    grid_->sendMessages(rho); //- For correct gradient computation

    //- Update the gravitational source term
    gradRho.computeFaces();
    sg.savePreviousTimeStep(timeStep, 1.);
    for (const Face &face: grid().faces())
        sg(face) = dot(g_, -face.centroid()) * gradRho(face);

    sg.oldField(0).faceToCell(rho, rho.oldField(0), fluid_);
    sg.faceToCell(rho, rho, fluid_);

    //- Update viscosity from kinematic viscosity
    mu.savePreviousTimeStep(timeStep, 1);
    mu.computeCells([this](const Cell &cell) {
        Scalar g = gamma(cell);
        return rho(cell) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
    });
    mu.computeFaces([this](const Face &face) {
        Scalar g = gamma(face);
        return rho(face) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
    });

    //- Update the surface tension
    ft.savePreviousTimeStep(timeStep, 1);
    ft.computeFaces(ib_);
    ft.oldField(0).faceToCell(rho, rho.oldField(0), fluid_);
    ft.faceToCell(rho, rho, fluid_);
}
