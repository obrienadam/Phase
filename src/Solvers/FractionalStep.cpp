#include "FractionalStep.h"
#include "FaceInterpolation.h"
#include "Source.h"
#include "SeoMittal.h"
#include "QuadraticIbm.h"

FractionalStep::FractionalStep(const Input &input,
                               std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        Solver(input, grid),
        u(addVectorField(input, "u")),
        p(addScalarField(input, "p")),
        gradP(addVectorField(std::make_shared<ScalarGradient>(p))),
        gradU(addTensorField(std::make_shared<JacobianField>(u))),
        uEqn_(input, u, "uEqn"),
        pEqn_(input, p, "pEqn"),
        fluid_(grid->createCellZone("fluid"))
{
    rho_ = input.caseInput().get<Scalar>("Properties.rho", 1);
    mu_ = input.caseInput().get<Scalar>("Properties.mu", 1);

    //- All active cells to fluid cells
    fluid_.add(grid_->localActiveCells());

    //- Create ib zones if any. Will also update local/global indices
    ib_.initCellZones(fluid_);
}

void FractionalStep::initialize()
{
    u.setBoundaryFaces();
    p.setBoundaryFaces();
}

std::string FractionalStep::info() const
{
    return "Fractional-step\n"
            "A simple 1-step fractional-step projection method\n"
            "May not produce accurate results near boundaries\n";
}

Scalar FractionalStep::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    //ib_.computeForce(rho_, mu_, u, p);
    //ib_.update(timeStep);

    printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    return 0;
}

Scalar FractionalStep::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for(const Cell& cell: fluid_)
    {
        Scalar co = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            co += std::max(dot(u(nb.face()), nb.outwardNorm()), 0.);

        for(const BoundaryLink& bd: cell.boundaries())
            co += std::max(dot(u(bd.face()), bd.outwardNorm()), 0.);

        co *= timeStep / cell.volume();
        maxCo = std::max(co, maxCo);
    }

    return grid_->comm().max(maxCo);
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

Scalar FractionalStep::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    //gradU.compute(fluid_);
    //grid_->sendMessages(gradU);

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.) + ib_.solidVelocity(u) == fv::laplacian(mu_/rho_, u, 0.5));
    //uEqn_ = (fv::ddt(u, timeStep) + qibm::div(u, u, ib_) + ib_.solidVelocity(u) == qibm::laplacian(mu_/rho_, u, ib_));
    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    u.interpolateFaces();

    //qibm::computeFaceVelocities(u, ib_);

    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    //pEqn_ = (seo::laplacian(ib_, rho_, timeStep, p) == seo::div(ib_, u));
    pEqn_ = (fv::laplacian(timeStep/rho_, p) + ib_.bcs(p) == src::div(u));
    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    //- Gradient
    p.setBoundaryFaces();
    gradP.computeFaces();

    return error;
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    //seo::correct(ib_, rho_, p, gradP, u, timeStep);
    //return;

    for(const Cell& cell: fluid_)
        u(cell) -= timeStep/rho_*gradP(cell);

    for(const Face& face: grid_->interiorFaces())
        u(face) -= timeStep/rho_*gradP(face);

    for(const Patch& patch: grid_->patches())
        switch(u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for(const Face& face: patch)
                    u(face) -= timeStep/rho_*gradP(face);
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for(const Face& face: patch)
                    u(face) = u(face.lCell()) - dot(u(face.lCell()), face.norm()) *face.norm()/face.norm().magSqr();
                break;
        }

    grid_->sendMessages(u);
}

Scalar FractionalStep::maxDivergenceError()
{
    Scalar maxError = 0.;

    for (const Cell &cell: fluid_)
    {
        Scalar div = 0.;

        for (const InteriorLink &nb: cell.neighbours())
            div += dot(u(nb.face()), nb.outwardNorm());

        for (const BoundaryLink &bd: cell.boundaries())
            div += dot(u(bd.face()), bd.outwardNorm());

        maxError = fabs(div) > maxError ? div : maxError;
    }

    return grid_->comm().max(maxError);
}


