#include "FractionalStepIncremental.h"
#include "ScalarGradient.h"
#include "Source.h"

FractionalStepIncremental::FractionalStepIncremental(const Input &input,
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
    rho_ = input.caseInput().get<Scalar>("Properties.rho", 1.);
    mu_ = input.caseInput().get<Scalar>("Properties.mu", 1.);

    //- All active cells to fluid cells
    fluid_.add(grid_->localActiveCells());

    //- Create ib zones if any. Will also update local/global indices
    ib_->initCellZones(fluid_);
}

void FractionalStepIncremental::initialize()
{
    u.setBoundaryFaces();
    p.setBoundaryFaces();
}

std::string FractionalStepIncremental::info() const
{
    return "Fractional-step\n"
            "An incremental fractional-step projection method\n";
}

Scalar FractionalStepIncremental::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    //ib_->update(timeStep);

    printf("Max Co = %lf\n", maxCourantNumber(timeStep));
    printf("Max divergence error = %.4e\n", maxDivergenceError());

    return 0.;
}

//- Protected methods

void FractionalStepIncremental::restartSolution()
{
    Solver::restartSolution();
    u.interpolateFaces();
    p.interpolateFaces();
}

Scalar FractionalStepIncremental::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(u, timeStep) + fv::div(u, u, 0.5) + ib_->velocityBcs(u)
             == fv::laplacian(mu_ / rho_, u, 0.5) - src::src(gradP / rho_, fluid_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    for (const Face &f: grid_->interiorFaces())
    {
        const Cell &l = f.lCell();
        const Cell &r = f.rCell();
        Scalar g = f.volumeWeight();

        u(f) = g * (u(l) + timeStep / rho_ * gradP(l))
               + (1. - g) * (u(r) + timeStep / rho_ * gradP(r))
               - timeStep / rho_ * gradP(f);
    }

    for (const Patch &patch: grid_->patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &face: patch)
                    u(face) = u(face.lCell()) + timeStep / rho_ * (gradP(face.lCell()) - gradP(face));
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &face: patch)
                    u(face) = u(face.lCell()) - dot(u(face.lCell()), face.norm()) * face.norm() / face.norm().magSqr();
                break;
        }

    return error;
}

Scalar FractionalStepIncremental::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p) + ib_->bcs(p) == src::div(u) + src::laplacian(timeStep / rho_, p));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);
    p.setBoundaryFaces();

    gradP.savePreviousTimeStep(timeStep, 1);
    gradP.compute(fluid_);
    grid_->sendMessages(gradP);

    return error;
}

void FractionalStepIncremental::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField &gradP0 = gradP.oldField(0);

    for (const Cell &cell: fluid_)
        u(cell) -= timeStep / rho_ * (gradP(cell) - gradP0(cell));

    grid_->sendMessages(u);

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho_ * (gradP(face) - gradP0(face));

    for (const Patch &patch: grid_->patches())
        switch (u.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &face: patch)
                    u(face) -= timeStep / rho_ * (gradP(face) - gradP0(face));
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &face: patch)
                    u(face) = u(face.lCell()) - dot(u(face.lCell()), face.norm()) * face.norm() / face.norm().magSqr();
                break;
        }
}

Scalar FractionalStepIncremental::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for (const Cell &cell: fluid_)
    {
        Scalar co = 0.;

        for (const InteriorLink &nb: cell.neighbours())
            co += std::max(dot(u(nb.face()), nb.outwardNorm()), 0.);

        for (const BoundaryLink &bd: cell.boundaries())
            co += std::max(dot(u(bd.face()), bd.outwardNorm()), 0.);

        maxCo = std::max(co * timeStep / cell.volume(), maxCo);
    }

    return grid_->comm().max(maxCo);
}

Scalar FractionalStepIncremental::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    Scalar co = maxCourantNumber(prevTimeStep);
    Scalar lambda1 = 0.1, lambda2 = 1.2;

    return grid_->comm().min(
            std::min(
                    std::min(maxCo / co * prevTimeStep, (1 + lambda1 * maxCo / co) * prevTimeStep),
                    std::min(lambda2 * prevTimeStep, maxTimeStep_)
            ));
}

Scalar FractionalStepIncremental::maxDivergenceError() const
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
