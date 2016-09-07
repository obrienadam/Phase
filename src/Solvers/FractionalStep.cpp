#include "FractionalStep.h"
#include "TotalVariationDiminishing.h"
#include "CrankNicolson.h"
#include "GradientEvaluation.h"
#include "FaceInterpolation.h"

FractionalStep::FractionalStep(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(addVectorField(input, "u")),
      sg(addVectorField("sg")),
      gradP(addVectorField("gradP")),
      p(addScalarField(input, "p")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      divUStar(addScalarField("uStar")),
      uEqn_(u, "uEqn", SparseMatrix::IncompleteLUT),
      pEqn_(p, "pEqn", SparseMatrix::IncompleteLUT)
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g"));

    volumeIntegrators_ = VolumeIntegrator::initVolumeIntegrators(input, *this);
    forceIntegrators_ = ForceIntegrator::initForceIntegrators(input, p, rho, mu, u);

    p.initNodes();

    uEqn_.matrix().setFill(1);
    pEqn_.matrix().setFill(8);

    interpolateFaces(fv::INVERSE_VOLUME, u);

    pEqn_ = (fv::laplacian(p) + ib_.eqns(p) == 0.);
    pEqn_.solve(); // dummy solve
}

std::string FractionalStep::info() const
{
    return Solver::info()
            + "Type: 2nd order fractional-step\n"
            + "Advection time-marching: Crank-Nicolson\n"
            + "Diffusion time-marching: Adams-Bashforth\n"
            + "Gravity: " + g_.toString() + "\n";
}

Scalar FractionalStep::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePEqn(timeStep);
    correctVelocity(timeStep);

    printf("Max Co = %lf\n", courantNumber(timeStep));

    return 0.;
}

Scalar FractionalStep::computeMaxTimeStep(Scalar maxCo) const
{
    Scalar maxTimeStepSqr = std::numeric_limits<Scalar>::infinity(), maxCoSqr = maxCo*maxCo;

    for(const Cell &cell: u.grid.fluidCells())
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaXSqr = nb.rCellVec().magSqr();
            Scalar magUSqr = u(nb.face()).magSqr();

            maxTimeStepSqr = std::min(maxTimeStepSqr, maxCoSqr*deltaXSqr/magUSqr);
        }

    return sqrt(maxTimeStepSqr);
}

//- Protected methods

Scalar FractionalStep::solveUEqn(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    sg.savePreviousTimeStep(timeStep, 1);
    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u, 1.5) + ib_.eqns(u) == cn::laplacian(mu, u, 0.5) - fv::source(gradP - sg));

    Scalar error = uEqn_.solve();

    computeFaceVelocities(timeStep);

    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    computeMassSource(timeStep);

    pEqn_ = (fv::laplacian(p) + ib_.eqns(p) == divUStar + fv::hydroStaticPressureBoundaries(p, rho, g_));

    Scalar error = pEqn_.solve(p.sparseVector(), false);

    for(const Face& face: p.grid.boundaryFaces())
    {
        if(p.boundaryType(face) == ScalarFiniteVolumeField::NORMAL_GRADIENT)
            p(face) = p(face.lCell()) + rho(face.lCell())*dot(g_, face.centroid() - face.lCell().centroid());
    }

    gradP.savePreviousTimeStep(timeStep, 1);
    computeGradient(fv::FACE_TO_CELL, p, gradP, true);

    return error;
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    const VectorFiniteVolumeField& gradP0 = gradP.prev(0);
    const VectorFiniteVolumeField& sg0 = sg.prev(0);

    for(const Cell &cell: grid_.fluidCells())
        u(cell) += timeStep/rho(cell)*(gradP0(cell) - sg0(cell) - gradP(cell) + sg(cell));

    for(const Face &face: grid_.interiorFaces()) // In the Francois paper, gradP0 is not part of the face correction
        u(face) += timeStep/rho(face)*(gradP0(face) - sg0(face) - gradP(face) + sg(face));

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
            u(face) += timeStep/rho(face)*(gradP0(face) - sg0(face) - gradP(face) + sg(face));
            break;
        };
    }
}

void FractionalStep::computeFaceVelocities(Scalar timeStep)
{
    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();
        const Scalar g = rCell.volume()/(rCell.volume() + lCell.volume());

        u(face) = g*(u(lCell) + timeStep/rho(lCell)*(gradP(lCell) - sg(lCell)))
                + (1. - g)*(u(rCell) + timeStep/rho(rCell)*(gradP(rCell) - sg(rCell)))
                - timeStep/rho(face)*(gradP(face) - sg(face));
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Cell &cell = face.lCell();

        switch(u.boundaryType(face))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) = u(cell) + timeStep/rho(cell)*(gradP(cell) - sg(cell))
                    - timeStep/rho(face)*(gradP(face) - sg(face));
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(cell.centroid());
            u(face) = u(cell) - dot(u(cell), nWall)*nWall/nWall.magSqr();
        }
            break;

        default:
            throw Exception("FractionalStep", "computeFaceVelocities", "unrecognized boundary condition type.");
        };
    }
}

void FractionalStep::computeMassSource(Scalar timeStep)
{
    divUStar.fill(0.);

    for(const Cell &cell: grid_.fluidCells())
    {
        for(const InteriorLink &nb: cell.neighbours())
            divUStar(cell) += rho(nb.face())/timeStep*dot(u(nb.face()), nb.outwardNorm()) + dot(gradP(nb.face()) - sg(nb.face()), nb.outwardNorm());

        for(const BoundaryLink &bd: cell.boundaries())
            divUStar(cell) += rho(bd.face())/timeStep*dot(u(bd.face()), bd.outwardNorm()) + dot(gradP(bd.face()) - sg(bd.face()), bd.outwardNorm());
    }
}

Scalar FractionalStep::courantNumber(Scalar timeStep)
{
    Scalar maxCoSqr = 0., timeStepSqr = timeStep*timeStep;

    for(const Cell &cell: u.grid.fluidCells())
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaXSqr = nb.rCellVec().magSqr();
            Scalar magUSqr = u(nb.face()).magSqr();

            maxCoSqr = std::max(maxCoSqr, magUSqr*timeStepSqr/deltaXSqr);
        }

    return sqrt(maxCoSqr);
}

