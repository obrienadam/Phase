#include "FractionalStep.h"
#include "CrankNicolson.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"
#include "SourceEvaluation.h"
#include "Source.h"
#include "GhostCellImmersedBoundary.h"
#include "MovingGhostCellImmersedBoundary.h"

FractionalStep::FractionalStep(const Input &input, const Communicator &comm, FiniteVolumeGrid2D &grid)
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
    alphaAdv_ = input.caseInput().get<Scalar>("Solver.CrankNicholsonAdvection", 0.5);
    alphaDiff_ = input.caseInput().get<Scalar>("Solver.CrankNicholsonDiffusion", 0.5);

    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));

    rho.savePreviousTimeStep(0., 1);
    mu.savePreviousTimeStep(0., 1);

    //- All active cells to fluid cells
    grid_.createCellZone("fluid", grid_.getCellIds(grid_.localActiveCells()));

    //- Create ib zones if any. Will also update local/global indices
    ibObjManager_.initCellZones();
}

void FractionalStep::initialize()
{
    //- Ensure computations start with a valid pressure field
    solvePEqn(1.);
    computeFaceVelocities(1.);
    solvePEqn(1.);

    //- Ensure the computation starts with a valid velocity field
    correctVelocity(1.);
}

std::string FractionalStep::info() const
{
    return Solver::info()
           + "Type: 2nd order fractional-step\n"
           + "Advection time-marching implicit weight: " + std::to_string(alphaAdv_) + "\n"
           + "Diffusion time-marching implicit weight: " + std::to_string(alphaDiff_) + "\n";
}

Scalar FractionalStep::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);
    //ib::correctVelocity(ibObjManager_, rho, p, gradP, u, timeStep);

    ibObjManager_.computeForce(rho, mu, u, p);
    ibObjManager_.update(timeStep);

    comm_.printf("Max Co = %lf\n", maxCourantNumber(timeStep));
    comm_.printf("Max divergence error = %.4e\n", maxDivergenceError());

    return 0.;
}

Scalar FractionalStep::maxCourantNumber(Scalar timeStep) const
{
    Scalar maxCo = 0;

    for (const Face &face: grid_.interiorFaces())
    {
        Vector2D sf = face.outwardNorm(face.lCell().centroid());
        Vector2D rc = face.rCell().centroid() - face.lCell().centroid();

        maxCo = std::max(maxCo, fabs(dot(u(face), sf) / dot(rc, sf)));
    }

    return comm_.max(maxCo * timeStep);
}

Scalar FractionalStep::computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const
{
    Scalar co = maxCourantNumber(prevTimeStep);
    Scalar lambda1 = 0.1, lambda2 = 1.2;

    return comm_.min(
            std::min(
                    std::min(maxCo / co * prevTimeStep, (1 + lambda1 * maxCo / co) * prevTimeStep),
                    std::min(lambda2 * prevTimeStep, maxTimeStep_)
            ));
}

//- Protected methods

Scalar FractionalStep::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);
    //uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho, u, u) + ib::gc(ibObjs(), u) ==
    //         fv::laplacian(mu, u));

    uEqn_ = ib::momentumEqn(ibObjManager_, rho, mu, p, u, timeStep);
    Scalar error = uEqn_.solve();
    grid_.sendMessages(comm_, u);
    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    pEqn_ = (fv::laplacian(timeStep / rho, p) + ib::gc(ibObjs(), p) ==
             source::div(u));

    //pEqn_ = ib::pressureEqn(ibObjManager_, rho, u, p, timeStep);
    Scalar error = pEqn_.solve();
    grid_.sendMessages(comm_, p);

    //- Compute pressure gradient
    gradP.savePreviousTimeStep(timeStep, 1);
    fv::computeGradient(fv::FACE_TO_CELL, p, gradP, false);
    grid_.sendMessages(comm_, gradP);

    return error;
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    for (const Cell &cell: grid_.cellZone("fluid"))
        u(cell) -= timeStep / rho(cell) * gradP(cell);

    grid_.sendMessages(comm_, u);

    for (const Face &face: grid_.interiorFaces())
        u(face) -= timeStep / rho(face) * gradP(face);

    for (const Face &face: grid_.boundaryFaces())
    {
        switch (u.boundaryType(face))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = face.outwardNorm(face.lCell().centroid());
                u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall) * nWall / nWall.magSqr();
            }
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                u(face) -= timeStep / rho(face) * gradP(face);
                break;
        };
    }
}

void FractionalStep::computeFaceVelocities(Scalar timeStep)
{
    for (const Face &face: u.grid.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();
        const Scalar g = rCell.volume() / (rCell.volume() + lCell.volume());

        u(face) = g * u(lCell) + (1. - g) * u(rCell);
    }

    for (const Face &face: u.grid.boundaryFaces())
    {
        const Cell &cell = face.lCell();

        switch (u.boundaryType(face))
        {
            case VectorFiniteVolumeField::FIXED:
                break;

            case VectorFiniteVolumeField::NORMAL_GRADIENT:
                u(face) = u(cell);
                break;

            case VectorFiniteVolumeField::SYMMETRY:
            {
                const Vector2D nWall = face.outwardNorm(cell.centroid());
                u(face) = u(cell) - dot(u(cell), nWall) * nWall / nWall.magSqr();
            }
                break;

                //default:
                //throw Exception("FractionalStep", "computeFaceVelocities", "unrecognized boundary condition type.");
        };
    }
}

Vector2D FractionalStep::maxVelocity() const
{
    Vector2D maxVelocity(0., 0.);
    for (const Cell &cell: grid_.cells())
        if (u(cell).magSqr() > maxVelocity.magSqr())
            maxVelocity = u(cell);

    return maxVelocity;
}

Vector2D FractionalStep::maxFaceVelocity() const
{
    Vector2D maxVelocity(0., 0.);
    for (const Face &face: grid_.faces())
        if (u(face).magSqr() > maxVelocity.magSqr())
            maxVelocity = u(face);

    return maxVelocity;
}

Scalar FractionalStep::maxDivergenceError() const
{
    Scalar maxError = 0.;

    for (const Cell &cell: grid_.cellZone("fluid"))
    {
        Scalar div = 0.;

        for (const InteriorLink &nb: cell.neighbours())
            div += dot(u(nb.face()), nb.outwardNorm());

        for (const BoundaryLink &bd: cell.boundaries())
            div += dot(u(bd.face()), bd.outwardNorm());

        if (fabs(div) > maxError)
            maxError = fabs(div);
    }

    return comm_.max(maxError);
}

