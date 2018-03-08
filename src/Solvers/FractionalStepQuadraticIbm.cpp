#include "FractionalStepQuadraticIbm.h"
#include "QuadraticIbm.h"
#include "Source.h"

FractionalStepQuadraticIbm::FractionalStepQuadraticIbm(const Input &input)
        :
        FractionalStep(input)
{
    for (const auto &ibObj: *ib_)
    {
        if (ibObj->type() != ImmersedBoundaryObject::QUADRATIC)
            throw Exception("FractionalStepQuadraticIbm",
                            "FractionalStepQuadraticIbm",
                            "immersed boundary object \"" + ibObj->name() + "\" is not type \"quadratic\".");
    }

    gradP.savePreviousTimeStep(0., 1);
}

Scalar FractionalStepQuadraticIbm::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(u, timeStep) + qibm::div(u, u, *ib_, 1.) + ib_->velocityBcs(u)
             == qibm::laplacian(mu_ / rho_, u, *ib_, 1.) - src::src(gradP / rho_, fluid_));

    Scalar error = uEqn_.solve();

    for(const Cell &cell: fluid_)
        u(cell) += timeStep / rho_ * gradP(cell);

    grid_->sendMessages(u);
    u.interpolateFaces();

//    for (const Face &f: grid_->interiorFaces())
//    {
//        Scalar g = f.volumeWeight();
//        const Cell &l = f.lCell();
//        const Cell &r = f.rCell();
//
//        u(f) = g * (u(l) - timeStep / rho_ * gradP0(l))
//               + (1. - g) * (u(r) - timeStep / rho_ * gradP0(r));
//    }
//
//    for (const Patch &patch: u.grid()->patches())
//        switch (u.boundaryType(patch))
//        {
//            case VectorFiniteVolumeField::FIXED:break;
//            case VectorFiniteVolumeField::NORMAL_GRADIENT:
//                for (const Face &f: patch)
//                {
//                    const Cell &l = f.lCell();
//                    u(f) = u(l) - timeStep / rho_ * gradP0(l);
//                }
//                break;
//            case VectorFiniteVolumeField::SYMMETRY:
//                for (const Face &f: patch)
//                    u(f) = u(f.lCell()) - dot(u(f.lCell()), f.norm()) * f.norm() / f.norm().magSqr();
//                break;
//        }

    return error;
}

Scalar FractionalStepQuadraticIbm::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho_, p, grid_->localActiveCells()) == src::div(u, grid_->localActiveCells()));

    Scalar error = pEqn_.solve();
    grid_->sendMessages(p);

    //- Gradient
    p.setBoundaryFaces();
    gradP.savePreviousTimeStep(timeStep, 1);
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