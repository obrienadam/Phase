#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"

#include "FractionalStep.h"

FractionalStep::FractionalStep(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
        :
        Solver(input, grid),
        fluid_(*cells_),
        u(*addField<Vector2D>(input, "u")),
        p(*addField<Scalar>(input, "p")),
        gradP(*std::static_pointer_cast<ScalarGradient>(addField<Vector2D>(std::make_shared<ScalarGradient>(p)))),
        gradU(*std::static_pointer_cast<JacobianField>(addField<Tensor2D>(std::make_shared<JacobianField>(u)))),
        uEqn_(input, u, "uEqn"),
        pEqn_(input, p, "pEqn")
{
    rho_ = input.caseInput().get<Scalar>("Properties.rho", 1);
    mu_ = input.caseInput().get<Scalar>("Properties.mu", 1);
    g_ = input.caseInput().get<std::string>("Properties.g", "(0,0)");

    //- All active cells to fluid cells
    //fluid_.add(grid_->localCells());
}

void FractionalStep::initialize()
{
    u.interpolateFaces();
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
    grid_->comm().printf("Updating IB positions...\n");
    ib_->update(timeStep);

    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    grid_->comm().printf("Max divergence error = %.4e\n", grid_->comm().max(maxDivergenceError()));
    grid_->comm().printf("Max CFL number = %.4lf\n", maxCourantNumber(timeStep));

    grid_->comm().printf("Computing IB forces...\n");
    ib_->computeForce(rho_, mu_, u, p, g_);

    return 0;
}

Scalar FractionalStep::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for (const Cell &cell: fluid_)
    {
        Scalar co = 0.;

        for (const InteriorLink &nb: cell.neighbours())
            co += std::max(dot(u(nb.face()), nb.outwardNorm()), 0.);

        for (const BoundaryLink &bd: cell.boundaries())
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

    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.) + ib_->velocityBcs(u)
             == fv::laplacian(mu_ / rho_, u, 0.5) - src::src(gradP / rho_));

    Scalar error = uEqn_.solve();

    for (const Cell &cell: fluid_)
        u(cell) += timeStep / rho_ * gradP(cell);

    grid_->sendMessages(u);
    u.interpolateFaces();

    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p) + ib_->bcs(p) == src::div(u));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    //- Gradient
    p.setBoundaryFaces();
    gradP.compute(fluid_);

    return error;
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: fluid_)
        u(cell) -= timeStep / rho_ * gradP(cell);

    grid_->sendMessages(u); //- Necessary

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho_ * gradP(face);

    for (const FaceGroup &patch: grid_->patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &face: patch)
                    u(face) -= timeStep / rho_ * gradP(face);
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &face: patch)
                    u(face) = u(face.lCell()) - dot(u(face.lCell()), face.norm()) * face.norm() / face.norm().magSqr();
                break;
        }
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