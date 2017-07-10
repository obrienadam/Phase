#include "FractionalStep.h"
#include "ScalarGradient.h"
#include "FaceInterpolation.h"
#include "Source.h"

FractionalStep::FractionalStep(const Input &input,
                               std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        Solver(input, grid),
        u(addVectorField(input, "u")),
        phi(addScalarField("phi")),
        p(addScalarField(input, "p")),
        gradP(addVectorField(std::make_shared<ScalarGradient>(p))),
        uEqn_(input, u, "uEqn"),
        pEqn_(input, phi, "pEqn"),
        fluid_(grid->createCellZone("fluid"))
{
    alphaAdv_ = input.caseInput().get<Scalar>("Solver.CrankNicholsonAdvection", 0.5);
    alphaDiff_ = input.caseInput().get<Scalar>("Solver.CrankNicholsonDiffusion", 0.5);

    rho_ = input.caseInput().get<Scalar>("Properties.rho", 1.);
    mu_ = input.caseInput().get<Scalar>("Properties.mu", 1.);

    phi.copyBoundaryTypes(p);
    ib_.copyBoundaryTypes(p, phi);

    //- All active cells to fluid cells
    fluid_.add(grid_->localActiveCells());

    //- Create ib zones if any. Will also update local/global indices
    ib_.initCellZones(fluid_);
}

void FractionalStep::initialize()
{
    //- Ensure computations start with a valid pressure field
    //    computeFaceVelocities(0.);
    //    solvePEqn(1.);

    //- Ensure the computation starts with a valid velocity field
    //    correctVelocity(1.);
}

std::string FractionalStep::info() const
{
    return Solver::info()
           + "Type: 2nd order fractional-step\n"
           + "Advection time-marching implicit weight: " + std::to_string(alphaAdv_) + "\n"
           + "Diffusion time-marching implicit weight: " + std::to_string(alphaDiff_) + "\n";
}

Scalar FractionalStep::solve(Scalar timeStep)
{
    solveUEqn(timeStep);

    ib_.clearFreshCells(); //- Clear the fresh cells

    solvePEqn(timeStep);
    correctVelocity(timeStep);

    //- ibObjManager_.computeForce(rho, mu, u, p);
    ib_.update(timeStep);

    printf("Max Co = %lf\n", maxCourantNumber(timeStep));
    printf("Max divergence error = %.4e\n", maxDivergenceError());

    return 0.;
}

//- Protected methods

Scalar FractionalStep::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho_, u, timeStep) + fv::div(rho_*u, u) + ib_.bcs(u) ==
             fv::laplacian(mu_, u) - fv::source(gradP, fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, phi) + ib_.bcs(phi) == fv::src::div(u));

    Scalar error = pEqn_.solve();

    for (const Cell &cell: grid_->localActiveCells())
        p(cell) += phi(cell);

    grid_->sendMessages(p);

    //- Compute pressure gradient
    p.interpolateFaces([](const Face &f) {
        return f.rCell().volume()/(f.lCell().volume() + f.rCell().volume());
    });

    gradP.savePreviousTimeStep(timeStep, 1);
    gradP.compute(fluid_);

    grid_->sendMessages(gradP);

    return error;
}

void FractionalStep::computeFaceVelocities(Scalar timeStep)
{
    auto alpha = [](const Face &f) {
        return f.rCell().volume()/(f.lCell().volume() + f.rCell().volume());
    };

    for (const Face &face: u.grid().interiorFaces())
    {
        Scalar g = alpha(face);

        u(face) = g * (u(face.lCell()) + timeStep / rho_ * gradP(face.lCell()))
                  + (1. - g) * (u(face.rCell()) + timeStep / rho_ * gradP(face.rCell()))
                  - timeStep / rho_ * gradP(face);
    }

    u.setBoundaryFaces(VectorFiniteVolumeField::SYMMETRY, [this](const Face& face){
        Vector2D nWall = face.outwardNorm(face.lCell().centroid());
        return u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
    });

    u.setBoundaryFaces(VectorFiniteVolumeField::NORMAL_GRADIENT, [this, timeStep](const Face& face){
        return u(face) + timeStep / rho_ * (gradP(face.lCell()) - gradP(face));
    });
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField &gradP0 = gradP.oldField(0);

    for (const Cell &cell: fluid_)
        u(cell) -= timeStep / rho_ * (gradP(cell) - gradP0(cell));

    grid_->sendMessages(u);

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho_ * (gradP(face) - gradP0(face));

    u.setBoundaryFaces(VectorFiniteVolumeField::SYMMETRY, [this](const Face& face){
        Vector2D nWall = face.outwardNorm(face.lCell().centroid());
        return u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
    });

    u.setBoundaryFaces(VectorFiniteVolumeField::NORMAL_GRADIENT, [this, &gradP0, timeStep](const Face& face){
        return u(face) - timeStep / rho_ * (gradP(face) - gradP0(face));
    });
}

Scalar FractionalStep::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for (const Face &face: grid_->interiorFaces())
    {
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();

        maxCo = std::max(maxCo, fabs(dot(u(face), sf) / dot(rc, sf)));
    }

    return grid_->comm().max(maxCo * timeStep);
}

Scalar FractionalStep::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    Scalar co = maxCourantNumber(prevTimeStep);
    Scalar lambda1 = 0.1, lambda2 = 1.2;

    return grid_->comm().min(
            std::min(
                    std::min(maxCo / co * prevTimeStep, (1 + lambda1 * maxCo / co) * prevTimeStep),
                    std::min(lambda2 * prevTimeStep, maxTimeStep_)
            ));
}

Scalar FractionalStep::maxDivergenceError() const
{
    Scalar maxError = 0.;

    for (const Cell &cell: fluid_)
    {
        Scalar div = 0.;

        for (const InteriorLink &nb: cell.neighbours())
            div += dot(u(nb.face()), nb.outwardNorm());

        for (const BoundaryLink &bd: cell.boundaries())
            div += dot(u(bd.face()), bd.outwardNorm());

        if (fabs(div) > maxError)
            maxError = fabs(div);
    }

    return grid_->comm().max(maxError);
}
