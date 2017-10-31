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
    for (auto ibObj: ib_.ibObjPtrs())
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("FractionalStepMultiphaseQuadraticIbm",
                            "FractionalStepMultiphaseQuadraticIbm",
                            "immersed boundary object \"" + ibObj->name() + "\" is not type \"quadratic\".");
    }

    addScalarField("ps");
}

Scalar FractionalStepMultiphaseQuadraticIbm::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep);
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    ib_.computeForce(rho, mu, u, p, g_);
    //ib_.computeForce(rho1_, mu1_, u, p, g_);
    ib_.update(timeStep);

    ScalarFiniteVolumeField& ps = scalarField("ps");
    for(const Cell& cell: grid_->cells())
        ps(cell) = p(cell) + rho(cell) * dot(g_, cell.centroid());

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

Scalar FractionalStepMultiphaseQuadraticIbm::solveGammaEqn(Scalar timeStep)
{
    auto beta = cicsam::beta(u, gradGamma, gamma, timeStep, 0.5);

    //- Advect volume fractions
    gamma.savePreviousTimeStep(timeStep, 1);
    gammaEqn_ = (fv::ddt(gamma, timeStep, fluid_) + cicsam::div(u, beta, gamma, fluid_, 0.5)
                 + ft.contactLineBcs(ib_) == 0.);

    Scalar error = gammaEqn_.solve();

    for(const Cell& cell: grid().localActiveCells())
        gamma(cell) = clamp(gamma(cell), 0., 1.);

    grid_->sendMessages(gamma);
    gamma.interpolateFaces();

    //- Update the gradient
    gradGamma.compute(fluid_);
    grid_->sendMessages(gradGamma);

    //- Update all other properties
    computeMomentumFlux(beta, timeStep);
    updateProperties(timeStep);

    return error;
}

Scalar FractionalStepMultiphaseQuadraticIbm::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + qibm::div(rhoU, u, ib_, 0.5) + ib_.velocityBcs(u)
             == qibm::laplacian(mu, u, ib_) + src::src(ft + sg, fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);
    //u.interpolateFaces();
    qibm::computeFaceVelocities(u, ib_);

//    gradP.faceToCell(rho, rho.oldField(0), grid_->localActiveCells());
//
//    for (const Face &f: grid_->interiorFaces())
//    {
//        Scalar g = f.volumeWeight();
//        const Cell &l = f.lCell();
//        const Cell &r = f.rCell();
//
//        u(f) = g * (u(l) - timeStep / rho(l) * gradP(l))
//               + (1. - g) * (u(r) - timeStep / rho(r) * gradP(r))
//               + timeStep / rho(f) * gradP(f);
//    }
//
//    for (const Patch &patch: u.grid().patches())
//        switch (u.boundaryType(patch))
//        {
//            case VectorFiniteVolumeField::FIXED:
//                break;
//            case VectorFiniteVolumeField::NORMAL_GRADIENT:
//                for (const Face &f: patch)
//                {
//                    const Cell &l = f.lCell();
//                    u(f) = u(l) - timeStep / rho(l) * gradP(l)
//                           + timeStep / rho(f) * gradP(f);
//                }
//                break;
//            case VectorFiniteVolumeField::SYMMETRY:
//                for (const Face &f: patch)
//                    u(f) = u(f.lCell()) - dot(u(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
//                break;
//        }

    return error;
}

Scalar FractionalStepMultiphaseQuadraticIbm::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho, p, grid_->localActiveCells()) == src::div(u, grid().localActiveCells()));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    p.setBoundaryFaces();
    gradP.computeFaces();
    gradP.faceToCell(rho, rho, grid().localActiveCells());

    return error;
}

void FractionalStepMultiphaseQuadraticIbm::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: grid().localActiveCells())
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