#include "FiniteVolume/Equation/AxisymmetricTimeDerivative.h"
#include "FiniteVolume/Equation/AxisymmetricDivergence.h"
#include "FiniteVolume/Equation/AxisymmetricLaplacian.h"
#include "FiniteVolume/Equation/AxisymmetricSource.h"

#include "FractionalStepAxisymmetric.h"

FractionalStepAxisymmetric::FractionalStepAxisymmetric(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
        :
        FractionalStep(input, grid)
{

}

Scalar FractionalStepAxisymmetric::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (axi::ddt(u_, timeStep) + axi::div(u_, u_, 1)== axi::vectorLaplacian(mu_ / rho_, u_, 1));
    Scalar error = uEqn_.solve();

    u_.interpolateFaces();

    return error;
}

Scalar FractionalStepAxisymmetric::solvePEqn(Scalar timeStep)
{
    pEqn_ = (axi::laplacian(timeStep / rho_, p_, 1.) == axi::src::div(u_));
    Scalar error = pEqn_.solve();

    p_.interpolateFaces();
    gradP_.computeAxisymmetric(*fluid_);

    return error;
}

void FractionalStepAxisymmetric::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: grid_->localCells())
        u_(cell) -= timeStep / rho_ * gradP_(cell);

    grid_->sendMessages(u_);

    for (const Face &face: grid_->interiorFaces())
        u_(face) -= timeStep / rho_ * gradP_(face);

    for (const FaceGroup &patch: grid_->patches())
        switch (u_.boundaryType(patch))
        {
            case VectorFiniteVolumeField::FIXED:
                break;
            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                for (const Face &face: patch)
                    u_(face) -= timeStep / rho_ * gradP_(face);
                break;
            case VectorFiniteVolumeField::SYMMETRY:
                for (const Face &face: patch)
                    u_(face) = u_(face.lCell()) - dot(u_(face.lCell()), face.norm()) * face.norm() / face.norm().magSqr();
                break;
        }
}

Scalar FractionalStepAxisymmetric::maxDivergenceError()
{
    Scalar maxError = 0.;

    for (const Cell &cell: *fluid_)
    {
        Scalar divU = 0.;
        for (const InteriorLink &nb: cell.neighbours())
            divU += dot(u_(nb.face()), nb.face().polarOutwardNorm(cell.centroid()));

        for (const BoundaryLink &bd: cell.boundaries())
            divU += dot(u_(bd.face()), bd.face().polarOutwardNorm(cell.centroid()));

        maxError = std::abs(divU) > maxError ? std::abs(divU) : maxError;
    }

    return grid_->comm().max(maxError);
}
