#include "FractionalStepAxisymmetric.h"
#include "AxisymmetricTimeDerivative.h"
#include "AxisymmetricDivergence.h"
#include "AxisymmetricLaplacian.h"
#include "AxisymmetricSource.h"

FractionalStepAxisymmetric::FractionalStepAxisymmetric(const Input &input,
                                                       std::shared_ptr<FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid)
{

}

Scalar FractionalStepAxisymmetric::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (axi::ddt(u, timeStep) + axi::div(u, u, 1)
             == axi::vectorLaplacian(mu_ / rho_, u, 1));
    Scalar error = uEqn_.solve();

    u.interpolateFaces();

    return error;
}

Scalar FractionalStepAxisymmetric::solvePEqn(Scalar timeStep)
{
    pEqn_ = (axi::laplacian(timeStep / rho_, p, 1.) == axi::src::div(u));
    Scalar error = pEqn_.solve();

    p.interpolateFaces();
    gradP.computeAxisymmetric(fluid_);

    return error;
}

void FractionalStepAxisymmetric::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: grid().localActiveCells())
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

Scalar FractionalStepAxisymmetric::maxDivergenceError()
{
    Scalar maxError = 0.;

    for (const Cell &cell: fluid_)
    {
        Scalar divU = 0.;
        for (const InteriorLink &nb: cell.neighbours())
            divU += dot(u(nb.face()), nb.face().polarOutwardNorm(cell.centroid()));

        for (const BoundaryLink &bd: cell.boundaries())
            divU += dot(u(bd.face()), bd.face().polarOutwardNorm(cell.centroid()));

        maxError = std::abs(divU) > maxError ? std::abs(divU) : maxError;
    }

    return grid_->comm().max(maxError);
}