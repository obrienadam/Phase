#include "FractionalStepExperimental.h"
#include "CrankNicolson.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"
#include "SourceEvaluation.h"
#include "Source.h"
#include "GhostCellImmersedBoundary.h"
//#include "MovingGhostCellImmersedBoundary.h"
#include "FiniteVolumeEquation.h"

FractionalStepExperimental::FractionalStepExperimental(const Input &input, const Communicator &comm, FiniteVolumeGrid2D& grid)
    :
      Solver(input, comm, grid),
      u(addVectorField(input, "u")),
      gradP(addVectorField("gradP")),
      p(addScalarField(input, "p")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      uEqn_(input, comm, u, "uEqn"),
      pEqn_(input, comm, p, "pEqn")
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));

    rho.savePreviousTimeStep(0., 1);
    mu.savePreviousTimeStep(0., 1);

    //- All active cells to fluid cells
    grid_.createCellZone("fluid", grid_.getCellIds(grid_.localActiveCells()));

    //- Create ib zones if any. Will also update local/global indices
    ibObjManager_.initCellZones();
}

std::string FractionalStepExperimental::info() const
{
    return Solver::info()
            + "Type: 2nd order fractional-step\n"
            + "Advection time-marching: Crank-Nicolson\n"
            + "Diffusion time-marching: Adams-Bashforth\n";
}

Scalar FractionalStepExperimental::solve(Scalar timeStep)
{
    ibObjManager_.update(timeStep);
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    comm_.printf("Max Co = %lf\n", maxCourantNumber(timeStep));

    return 0.;
}

Scalar FractionalStepExperimental::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for(const Face &face: grid_.interiorFaces())
    {
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();

        maxCo = std::max(maxCo, fabs(dot(u(face), sf)/dot(rc, sf)));
    }

    return comm_.max(maxCo*timeStep);
}

Scalar FractionalStepExperimental::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    Scalar co = maxCourantNumber(prevTimeStep);
    Scalar lambda1 = 0.1, lambda2 = 1.2;

    return comm_.min(
                std::min(
                    std::min(maxCo/co*prevTimeStep, (1 + lambda1*maxCo/co)*prevTimeStep),
                    std::min(lambda2*prevTimeStep, maxTimeStep_)
                    ));
}

//- Protected methods

Scalar FractionalStepExperimental::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho, u, u, 0.5) + ib::gc(ibObjs(), u) == cn::laplacian(mu, u, 0.5));

    Scalar error = uEqn_.solve();

    grid_.sendMessages(comm_, u);

    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStepExperimental::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep/rho, p) + ib::gc(ibObjs(), p) == source::div(u));
    Scalar error = pEqn_.solveWithGuess();

    grid_.sendMessages(comm_, p);

    p.setBoundaryFaces();
    fv::computeInverseWeightedGradient(rho, p, gradP);

    grid_.sendMessages(comm_, gradP);

    return error;
}

void FractionalStepExperimental::correctVelocity(Scalar timeStep)
{
    for(const Cell &cell: grid_.cellZone("fluid"))
        u(cell) -= timeStep/rho(cell)*gradP(cell);

    grid_.sendMessages(comm_, u);

    for(const Face &face: grid_.interiorFaces())
        u(face) -= timeStep/rho(face)*gradP(face);

    for(const Face &face: grid_.boundaryFaces())
    {
        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall)*nWall/nWall.magSqr();
        }
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) -= timeStep/rho(face)*gradP(face);
            break;
        };
    }
}

void FractionalStepExperimental::computeFaceVelocities(Scalar timeStep)
{
    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();
        const Scalar g = rCell.volume()/(rCell.volume() + lCell.volume());

        u(face) = g*u(lCell) + (1. - g)*u(rCell);
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Cell &cell = face.lCell();

        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) = u(cell);
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(cell.centroid());
            u(face) = u(cell) - dot(u(cell), nWall)*nWall/nWall.magSqr();
        }
            break;

            //default:
            //throw Exception("FractionalStepExperimental", "computeFaceVelocities", "unrecognized boundary condition type.");
        };
    }
}

