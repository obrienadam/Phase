#include "Algorithm.h"
#include "FractionalStepMultiphase.h"
#include "CrankNicolson.h"
#include "Cicsam.h"
#include "FaceInterpolation.h"
#include "Source.h"

FractionalStepMultiphase::FractionalStepMultiphase(const Input &input,
                                                   std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid),
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
    addVectorField(ft.nPtr());

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

void FractionalStepMultiphase::initialize()
{
    FractionalStep::initialize();

    //- Ensure the computation starts with a valid gamma field
    gradGamma.compute(fluid_);

    //- Call twice to init a time history
    gamma.savePreviousTimeStep(0, 1);
    updateProperties(0.);
    updateProperties(0.);
}

Scalar FractionalStepMultiphase::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep); //- Solve a sharp gamma equation
    solveUEqn(timeStep); //- Solve a momentum prediction
    solvePEqn(timeStep); //- Solve the pressure equation, using sharp value of rho
    correctVelocity(timeStep);

    printf("Max Co = %lf\n", maxCourantNumber(timeStep));
    printf("Max absolute velocity divergence error = %.4e\n", maxDivergenceError());

    return 0.;
}

Scalar FractionalStepMultiphase::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    return std::min( //- Both args have already been globalling minimized
            FractionalStep::computeMaxTimeStep(maxCo, prevTimeStep),
            capillaryTimeStep_
    );
}

Scalar FractionalStepMultiphase::solveGammaEqn(Scalar timeStep)
{
    gamma.savePreviousTimeStep(timeStep, 1);
    auto beta = cicsam::computeBeta(u, gradGamma, gamma, timeStep, cicsamBlending_);
    gammaEqn_ = (fv::ddt(gamma, timeStep) + cicsam::div(u, beta, gamma) + ib_.bcs(gamma) == 0.);
    Scalar error = gammaEqn_.solve();

    // - While this may affect mass conservation, it prevents issues at high density ratios
    gamma.computeCells([this](const Cell &cell){
        return clamp(gamma(cell), 0., 1.);
    });

    grid_->sendMessages(gamma);
    gamma.interpolateFaces([this, &beta](const Face& face){
        return dot(u(face), face.outwardNorm(face.lCell().centroid())) > 0 ? 1. - beta(face): beta(face);
    });

    //- Recompute gradGamma for next cicsam iteration
    gamma.setBoundaryFaces();
    gradGamma.compute(fluid_);

    grid_->sendMessages(gradGamma); //- In case donor cell is on another proc
    updateProperties(timeStep);

    return error;
}

//- Protected methods
Scalar FractionalStepMultiphase::solveUEqn(Scalar timeStep)
{
    const ScalarFiniteVolumeField &rho0 = rho.oldField(0);
    const VectorFiniteVolumeField &sg0 = sg.oldField(0);
    const VectorFiniteVolumeField &ft0 = ft.oldField(0);

    gradU.compute(fluid_);
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u) + ib_.solidVelocity(u)
             == cn::laplacian(mu, u, 0.5) - src::ftc(rho, rho0, gradP - ft0 - sg0, fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    for(const Face& f: grid_->interiorFaces())
    {
        Scalar g = f.volumeWeight();
        const Cell& l = f.lCell();
        const Cell& r = f.rCell();

        u(f) = g*(u(l) + timeStep/rho0(l)*(gradP(l) - sg0(l) - ft0(l)))
               + (1. - g)*(u(r) + timeStep/rho0(r)*(gradP(r) - sg0(r) - ft0(r)))
               - timeStep/rho(f)*(gradP(f) - sg(f) - ft(f));
    }

    for(const Patch& patch: grid_->patches())
        switch(u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for(const Face& f: patch)
                {
                    const Cell& l = f.lCell();
                    u(f) = u(l) + timeStep/rho0(l)*(gradP(l) - sg0(l) - ft0(l))
                           - timeStep/rho(f)*(gradP(f) - sg(f) - ft(f));
                }
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for(const Face& f: patch)
                {
                    Vector2D t = f.norm().tangentVec();
                    u(f) = u(f.lCell()) - dot(u(f.lCell()), t)*t/t.magSqr();
                }
                break;
        }

    return error;
}

Scalar FractionalStepMultiphase::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho, p) + ib_.bcs(p) == src::div(u) + src::laplacian(timeStep / rho, p));
    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    //- Compute pressure gradient. Weighted gradients greatly reduce the effect of large pressure differences
    p.setBoundaryFaces();
    gradP.savePreviousTimeStep(timeStep, 1);
    gradP.compute(fluid_);
    grid_->sendMessages(gradP);

    return error;
}

void FractionalStepMultiphase::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField &gradP0 = gradP.oldField(0);
    const VectorFiniteVolumeField &sg0 = sg.oldField(0);
    const VectorFiniteVolumeField &ft0 = ft.oldField(0);
    const ScalarFiniteVolumeField &rho0 = rho.oldField(0);

    auto fb = src::ftc(rho0, rho0, gradP0 - sg0 - ft0, fluid_) / rho0
              - src::ftc(rho, rho, gradP - sg - ft, fluid_) / rho;

    for(const Cell& cell: fluid_)
        u(cell) += timeStep * fb(cell) / cell.volume();

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho(face) * (gradP(face) - gradP0(face));

    for(const Patch& patch: grid_->patches())
        switch(u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for(const Face& face: patch)
                    u(face) -= timeStep/rho(face)*(gradP(face) - gradP0(face));
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for(const Face& face: patch)
                {
                    Vector2D t = face.norm().tangentVec();
                    u(face) = u(face.lCell()) - dot(u(face.lCell()), t)*t/t.magSqr();
                }
                break;
        }

    grid_->sendMessages(u);
}

void FractionalStepMultiphase::updateProperties(Scalar timeStep)
{
    auto alpha = [](const Face& f){
        return f.volumeWeight();
    };

    //- Update rhoU
    rhoU.computeFaces([this](const Face& face){
        Scalar g = clamp(gamma(face), 0., 1.);
        Scalar g0 = clamp(gamma.oldField(0)(face), 0., 1.);
        Scalar rhoF = (1. - g)*rho1_ + g*rho2_;
        Scalar rhoF0 = (1. - g0)*rho1_ + g0*rho2_;
        return (rhoF + rhoF0)/2.*u(face);
    });

    //- Update density
    rho.savePreviousTimeStep(timeStep, 1);
    rho.computeCells([this](const Cell& cell){
        Scalar g = clamp(gamma(cell), 0., 1.);
        return (1. - g)*rho1_ + g*rho2_;
    });
    rho.interpolateFaces(alpha);

    //- Update viscosity
    mu.savePreviousTimeStep(timeStep, 1);
    mu.computeCells([this](const Cell& cell){
        Scalar g = clamp(gamma(cell), 0., 1.);
        return (1. - g)*mu1_ + g*mu2_;
    });
    mu.interpolateFaces(alpha);

    //- Update the surface tension
    ft.savePreviousTimeStep(timeStep, 1);
    //ft.compute();

    //- Update the gravitational source term
    gradRho.compute(fluid_);

    sg.savePreviousTimeStep(timeStep, 1);
    for(const Cell& cell: fluid_)
        sg(cell) = dot(g_, -cell.centroid()) * gradRho(cell);

    for(const Face& face: grid().faces())
        sg(face) = dot(g_, -face.centroid()) * gradRho(face);
}
