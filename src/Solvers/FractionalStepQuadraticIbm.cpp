#include "FractionalStepQuadraticIbm.h"
#include "QuadraticIbm.h"
#include "Source.h"

FractionalStepQuadraticIbm::FractionalStepQuadraticIbm(const Input &input, std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid)
{
    for (const auto &ibObj: *ib_)
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("FractionalStepQuadraticIbm",
                            "FractionalStepQuadraticIbm",
                            "immersed boundary object \"" + ibObj->name() + "\" is not type \"quadratic\".");
    }
}

Scalar FractionalStepQuadraticIbm::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    //gradU.compute(fluid_);
    //grid_->sendMessages(gradU);

    uEqn_ = (fv::ddt(u, timeStep) + qibm::div(u, u, *ib_) + ib_->velocityBcs(u) == qibm::laplacian(mu_ / rho_, u, *ib_));

    Scalar error = uEqn_.solve();
    grid_->sendMessages(u);

    u.interpolateFaces();
    //qibm::computeFaceVelocities(u, ib_);

    return error;
}

Scalar FractionalStepQuadraticIbm::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p, grid_->localActiveCells()) == src::div(u, grid_->localActiveCells()));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    //- Gradient
    p.setBoundaryFaces();
    gradP.compute(grid_->localActiveCells());

    return error;
}

void FractionalStepQuadraticIbm::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: grid_->localActiveCells())
        u(cell) -= timeStep / rho_ * gradP(cell);

    grid_->sendMessages(u);

    for (const Face &face: grid_->interiorFaces())
        u(face) -= timeStep / rho_ * gradP(face);

    for (const Patch &patch: grid_->patches())
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