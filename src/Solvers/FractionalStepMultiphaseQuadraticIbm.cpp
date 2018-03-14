#include "Algorithm.h"
#include "FractionalStepMultiphaseQuadraticIbm.h"
#include "Cicsam.h"
#include "QuadraticIbm.h"
#include "TimeDerivative.h"
#include "Divergence.h"
#include "Laplacian.h"
#include "Source.h"

FractionalStepMultiphaseQuadraticIbm::FractionalStepMultiphaseQuadraticIbm(const Input &input)
        :
        FractionalStepMultiphase(input),
        gammaCl(addScalarField("gammaCl")),
        gradGammaCl(addVectorField(std::make_shared<ScalarGradient>(gammaCl)))
{
    for (const auto &ibObj: *ib_)
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("FractionalStepMultiphaseQuadraticIbm",
                            "FractionalStepMultiphaseQuadraticIbm",
                            "immersed boundary object \"" + ibObj->name() + "\" is not of type \"quadratic\".");
    }

    gammaCl.copyBoundaryTypes(gamma);
}

Scalar FractionalStepMultiphaseQuadraticIbm::solve(Scalar timeStep)
{
    solveGammaEqn(timeStep);
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    ib_->computeForce(rho, mu, u, p, gamma, ft,  g_);
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
    gammaEqn_ = (fv::ddt(gamma, timeStep, fluid_) + ft.contactLineBcs(gamma)
                 + cicsam::div(u, beta, gamma, fluid_, 1)
                 == 0);

    Scalar error = gammaEqn_.solve();

    std::for_each(gamma.begin(), gamma.end(), [](Scalar &g)
    {
        g = clamp(g, 0., 1.);
    });

    grid_->sendMessages(gamma);
    gamma.interpolateFaces();

    //- Update the gradient
    gradGamma.compute(fluid_);
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
             == qibm::laplacian(mu, u, *ib_, 1) + src::src(ft + sg, fluid_));

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
            case VectorFiniteVolumeField::FIXED:break;
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
            case VectorFiniteVolumeField::FIXED:break;
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
    /*
    Equation<Scalar> eqn(gammaCl);
    eqn.setSparseSolver(std::make_shared<TrilinosAmesosSparseMatrixSolver>(grid_->comm()));

    for (const Cell &cell: fluid_)
    {
        eqn.set(cell, cell, 1.);
        eqn.setSource(cell, -gamma(cell));
    }

    for (const auto &ibObj: *ib_)
    {
        Scalar theta = ft.theta(*ibObj);

        for (const Cell &cell: ibObj->ibCells())
        {
            Vector2D wn = -ibObj->nearestEdgeNormal(cell.centroid());

            Ray2D r1 = Ray2D(cell.centroid(), wn.rotate(M_PI_2 - theta));
            Ray2D r2 = Ray2D(cell.centroid(), wn.rotate(theta - M_PI_2));

            GhostCellStencil m1(cell, ibObj->nearestIntersect(cell.centroid()), r1.r(), *grid_);
            GhostCellStencil m2(cell, ibObj->nearestIntersect(cell.centroid()), r2.r(), *grid_);

            Vector2D grad1 = m1.bpGrad(gamma);
            Vector2D grad2 = m2.bpGrad(gamma);

            if (dot(grad2, r2.r()) < dot(grad1, r1.r()))
                std::swap(m1, m2);

            if (theta > M_PI_2)
                eqn.add(m1.cell(), m1.neumannCells(), m1.neumannCoeffs());
            else
                eqn.add(m2.cell(), m2.neumannCells(), m2.neumannCoeffs());
        }

        for (const Cell &cell: ibObj->solidCells())
        {
            eqn.set(cell, cell, 1.);
            eqn.setSource(cell, 0.);
        }
    }

    eqn.solve();

    std::for_each(gammaCl.begin(), gammaCl.end(), [](Scalar &g)
    {
        g = clamp(g, 0., 1.);
    });

    grid_->sendMessages(gammaCl);

    gammaCl.setBoundaryFaces();
    gradGammaCl.compute(fluid_);
*/
    //- Update density
    rho.savePreviousTimeStep(timeStep, 1);
    rho.computeCells([this](const Cell &cell)
                     {
                         Scalar g = gamma(cell);
                         return (1. - g) * rho1_ + g * rho2_;
                     });

    // grid_->sendMessages(rho); //- For correct gradient computation

    rho.computeFaces([this](const Face &face)
                     {
                         Scalar g = gamma(face);
                         return (1. - g) * rho1_ + g * rho2_;
                     });

    //- Update the gravitational source term
    gradRho.computeFaces();
    sg.savePreviousTimeStep(timeStep, 1.);
    for (const Face &face: grid_->faces())
        sg(face) = dot(g_, -face.centroid()) * gradRho(face);

    sg.oldField(0).faceToCell(rho, rho.oldField(0), fluid_);
    sg.faceToCell(rho, rho, fluid_);

    //- Must be communicated for proper momentum interpolation
    grid_->sendMessages(sg.oldField(0));
    grid_->sendMessages(sg);

    //- Update viscosity from kinematic viscosity
    mu.savePreviousTimeStep(timeStep, 1);
    mu.computeCells([this](const Cell &cell)
                    {
                        Scalar g = gamma(cell);
                        return rho(cell) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
                    });

    // grid_->sendMessages(mu);

    mu.computeFaces([this](const Face &face)
                    {
                        Scalar g = gamma(face);
                        return rho(face) / ((1. - g) * rho1_ / mu1_ + g * rho2_ / mu2_);
                    });

    //- Update the surface tension
    ft.savePreviousTimeStep(timeStep, 1);
    ft.computeFaceInterfaceForces(gamma, gradGamma);
//
//    //- Predicate ensures cell-centred values aren't overwritten for cells neighbouring ib cells
//    auto p = [this](const Cell &cell)
//    {
//        for (const CellLink &nb: cell.neighbours())
//            if (ib_->ibObj(nb.cell().centroid()))
//                return false;
//        return true;
//    };

    ft.oldField(0).faceToCell(rho, rho.oldField(0), fluid_);
    ft.faceToCell(rho, rho, fluid_);

    //- Must be communicated for proper momentum interpolation
    grid_->sendMessages(ft.oldField(0));
    grid_->sendMessages(ft);
}