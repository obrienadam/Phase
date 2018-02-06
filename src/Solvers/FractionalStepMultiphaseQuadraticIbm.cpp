#include <Algorithm.h>
#include "FractionalStepMultiphaseQuadraticIbm.h"
#include "Cicsam.h"
#include "QuadraticIbm.h"
#include "Source.h"

FractionalStepMultiphaseQuadraticIbm::FractionalStepMultiphaseQuadraticIbm(const Input &input,
                                                                           std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStepMultiphase(input, grid)
{
    for (const auto &ibObj: *ib_)
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("FractionalStepMultiphaseQuadraticIbm",
                            "FractionalStepMultiphaseQuadraticIbm",
                            "immersed boundary object \"" + ibObj->name() + "\" is not type \"quadratic\".");
    }
}

Scalar FractionalStepMultiphaseQuadraticIbm::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep);
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    ib_->computeForce(rho, mu, u, p, g_);
    ib_->update(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

Scalar FractionalStepMultiphaseQuadraticIbm::solveGammaEqn(Scalar timeStep)
{
    cicsam::beta(u, gradGamma, gamma, timeStep, beta, 0.5);

    //- Advect volume fractions
    gamma.savePreviousTimeStep(timeStep, 1);
    gammaEqn_ = (fv::ddt(gamma, timeStep, grid_->localActiveCells()) + cicsam::div(u, beta, gamma, grid_->localActiveCells(), 1)
                 == ft.contactLineBcs(gamma));

    Scalar error = gammaEqn_.solve();

    std::for_each(gamma.begin(), gamma.end(), [](Scalar &g)
    {
        g = std::max(std::min(g, 1.), 0.);
    });

    grid_->sendMessages(gamma);
    gamma.interpolateFaces();

    //- Update the gradient
    gradGamma.compute(grid_->localActiveCells(), ScalarGradient::FACE_TO_CELL);
    grid_->sendMessages(gradGamma);

    //- Update all other properties
    cicsam::computeMomentumFlux(rho1_, rho2_, u, gamma, beta, timeStep, rhoU);
    updateProperties(timeStep);

    return error;
}

Scalar FractionalStepMultiphaseQuadraticIbm::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + qibm::div(rhoU, u, *ib_, 1) + ib_->velocityBcs(u)
             == qibm::laplacian(mu, u, *ib_) + src::src(ft + sg, fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    for (const Face &f: grid_->interiorFaces())
    {
        Scalar g = f.volumeWeight();
        const Cell &l = f.lCell();
        const Cell &r = f.rCell();

        u(f) = g * (u(l) - timeStep / rho(l) * (ft(l) + sg(l)))
               + (1. - g) * (u(r) - timeStep / rho(r) * (ft(r) + sg(r)))
               + timeStep / rho(f) * (ft(f) + sg(f));
    }

    for (const Patch &patch: u.grid()->patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &f: patch)
                {
                    const Cell &l = f.lCell();
                    u(f) = u(l) - timeStep / rho(l) * (ft(l) + sg(l))
                           + timeStep / rho(f) * (ft(f) + sg(f));
                }
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &f: patch)
                    u(f) = u(f.lCell()) - dot(u(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
                break;
        }

    return error;
}

Scalar FractionalStepMultiphaseQuadraticIbm::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho, p, grid_->localActiveCells()) == src::div(u, grid_->localActiveCells()));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    p.setBoundaryFaces();
    gradP.computeFaces();
    gradP.faceToCell(rho, rho, grid_->localActiveCells());

    return error;
}

void FractionalStepMultiphaseQuadraticIbm::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: grid_->localActiveCells()) // Try correcting over the entire domain!!
        u(cell) -= timeStep / rho(cell) * gradP(cell);

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho(face) * gradP(face);

    for (const Patch &patch: grid_->patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &face: patch)
                    u(face) -= timeStep / rho(face) * gradP(face);
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &f: patch)
                    u(f) = u(f.lCell()) - dot(u(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
                break;
        }

    grid_->sendMessages(u);
}

void FractionalStepMultiphaseQuadraticIbm::updateProperties(Scalar timeStep)
{
    FractionalStepMultiphase::updateProperties(timeStep);
}