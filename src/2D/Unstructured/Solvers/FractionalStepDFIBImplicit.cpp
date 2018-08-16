#include "FiniteVolume/Equation/TimeDerivative.h"
#include "FiniteVolume/Equation/Divergence.h"
#include "FiniteVolume/Equation/Laplacian.h"
#include "FiniteVolume/Equation/Source.h"

#include "FractionalStepDFIBImplicit.h"

FractionalStepDFIBImplicit::FractionalStepDFIBImplicit(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
    :
      FractionalStepDFIB(input, grid)
{

}

Scalar FractionalStepDFIBImplicit::solveUEqn(Scalar timeStep)
{
    u_.savePreviousTimeStep(timeStep, 1);

    uEqn_ = (fv::ddt(u_, timeStep) + fv::div(u_, u_, 0.)
             == fv::laplacian(mu_ / rho_, u_, 0.) - src::src(gradP_ / rho_) + ib_->velocityBcs(u_, timeStep));

    Scalar error = uEqn_.solve();

    Vector2D fib(0., 0.);

    for(const Cell& cell: *fluid_)
    {
        Vector2D f(0., 0.);

        f += rho_ * (u_(cell) - u_.oldField(0)(cell)) * cell.volume() / timeStep;

        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar flux = rho_ * dot(u_(nb.face()), nb.outwardNorm());
            f += std::max(flux, 0.) * u_.oldField(0)(cell);
            f += std::min(flux, 0.) * u_.oldField(0)(nb.cell());

            f -= mu_ * (u_.oldField(0)(nb.cell()) - u_.oldField(0)(cell)) * dot(nb.outwardNorm(), nb.rCellVec()) / nb.rCellVec().magSqr();
        }

        f += gradP_(cell) * cell.volume();


        fib += f;
    }

    fib = Vector2D(grid_->comm().sum(fib.x), grid_->comm().sum(fib.y));

    if(grid_->comm().isMainProc())
        std::cout << "FIB = " << fib << "\n";

    for(const Cell &c: u_.cells())
        u_(c) += timeStep / rho_ * gradP_(c);

    grid_->sendMessages(u_);
    u_.interpolateFaces();

//    const auto &u0 = u_.oldField(0);
//    for(const Face &f: grid_->interiorFaces())
//    {
//        const Cell &l = f.lCell();
//        const Cell &r = f.rCell();
//        Scalar g = f.volumeWeight();

//        u_(f) = g * (u_(l) - u0(l)) + (1. - g) * (u_(r) - u0(r)) + u0(f);
//    }

//    u_.setBoundaryFaces();

    return error;
}
