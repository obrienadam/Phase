#include "FractionalStepDirectForcingMultiphase.h"
#include "FiniteVolume/ImmersedBoundary/DirectForcingImmersedBoundaryObject.h"
#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"

FractionalStepDirectForcingMultiphase::FractionalStepDirectForcingMultiphase(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStepMultiphase(input, grid),
      fb(*addField<Vector2D>("fb"))
{

}

Scalar FractionalStepDirectForcingMultiphase::solve(Scalar timeStep)
{

}

Scalar FractionalStepDirectForcingMultiphase::solveUEqn(Scalar timeStep)
{
    //- Predictor step
    u.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u, 0.5) + ib_->velocityBcs(u)
             == fv::laplacian(mu, u, 0.5) + src::src(ft + sg - gradP));

    Scalar error = uEqn_.solve();

    fb.fill(Vector2D(0., 0.));

    for (const auto &ibObj: *ib_)
    {
        for(const Cell& cell: ibObj->ibCells())
        {
            DirectForcingImmersedBoundaryObject::Stencil st(u, cell, *ibObj);
            fb(cell) = (st.uf() - u(cell)) / timeStep;
        }

        for(const Cell& cell: ibObj->solidCells())
            fb(cell) = (ibObj->velocity(cell.centroid()) - u(cell)) / timeStep;
    }

    for (const Cell &cell: fluid_)
        u(cell) = u.oldField(0)(cell);

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rhoU, u, 0.5) + ib_->velocityBcs(u)
             == fv::laplacian(mu, u, 0.5) + src::src(ft + sg + fb - gradP));


    grid_->sendMessages(u);

    for (const Face &face: grid_->interiorFaces())
    {
        Scalar g = face.volumeWeight();
        const Cell& l = face.lCell();
        const Cell& r = face.rCell();

        u(face) = g * (u(l) - timeStep / rho(l) * (ft(l) + sg(l) - gradP(l)))
                + (1. - g) * (u(r) - timeStep / rho(r) * (ft(r) + sg(r) - gradP(r)))
                + timeStep / rho(face) * (ft(face) + sg(face));
    }

    for (const FaceGroup &patch: grid_->patches())
        switch (u.boundaryType(patch))
        {
        case VectorFiniteVolumeField::FIXED:
            break;
        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            for (const Face &f: patch)
            {
                const Cell &l = f.lCell();
                u(f) = u(l) - timeStep / rho(l) * (ft(l) + sg(l) - gradP(l))
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
