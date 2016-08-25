#include "FractionalStep.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"
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

    uEqn_.matrix().setFill(2);
    pEqn_.matrix().setFill(3);

    u.savePreviousTimeStep(0., 1);
    computeAdvectingVelocity(0.);
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

    for(const ForceIntegrator &fi: forceIntegrators_)
        fi.integrate();

    for(const VolumeIntegrator &vi: volumeIntegrators_)
        vi.integrate();

    printf("Current time step = %lf\n", timeStep);
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

    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u) + ib_.eqns(u) == ab::laplacian(mu, u) + fv::source(sg - gradP));

    Scalar error = uEqn_.solve();
    computeAdvectingVelocity(timeStep);

    return error;
}

Scalar FractionalStep::solvePEqn(Scalar timeStep)
{
    computeMassSource(timeStep);

    pEqn_ = (fv::laplacian(p) + ib_.eqns(p) == divUStar);
    Scalar error = pEqn_.solve(p.sparseVector());

    gradP.savePreviousTimeStep(timeStep, 1);

    computeGradient(fv::GREEN_GAUSS_CELL_CENTERED, p, gradP);

    return error;
}

void FractionalStep::correctVelocity(Scalar timeStep)
{
    const ScalarFiniteVolumeField &p0 = p.prev(0);
    const VectorFiniteVolumeField &gradP0 = gradP.prev(0);


    for(const Cell &cell: grid_.fluidCells())
        u(cell) += timeStep/rho(cell)*(gradP0(cell) - gradP(cell));

    for(const Face &face: grid_.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();
        Vector2D rc = rCell.centroid() - lCell.centroid();

        u(face) += timeStep/rho(face)*((p0(rCell) - p0(lCell) - p(rCell) + p(lCell))*rc/rc.magSqr());
    }

    for(const Face &face: grid_.boundaryFaces())
    {
        const Cell &cell = face.lCell();
        Vector2D rf = face.centroid() - cell.centroid();

        if(u.boundaryType(face.id()) == VectorFiniteVolumeField::FIXED)
            continue;
        else if(u.boundaryType(face.id()) == VectorFiniteVolumeField::NORMAL_GRADIENT) // Not correct
            continue;
        else if(u.boundaryType(face.id()) == VectorFiniteVolumeField::SYMMETRY) // This is not correct
            continue;
        else
            u(face) += timeStep/rho(face)*((p0(face) - p0(cell) - p(face) + p(cell))*rf/rf.magSqr());
    }
}

void FractionalStep::computeMassSource(Scalar timeStep)
{
    p.savePreviousTimeStep(timeStep, 1);
    divUStar.fill(0.);

    for(const Cell &cell: grid_.fluidCells())
    {
        for(const InteriorLink &nb: cell.neighbours())
        {
            const Vector2D& sf = nb.outwardNorm();
            const Vector2D& rc = nb.rCellVec();

            divUStar(cell) += rho(nb.face())/timeStep*dot(u(nb.face()), sf) + (p(nb.cell()) - p(cell))*dot(rc, sf)/dot(rc, rc);
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            const Vector2D& sf = bd.outwardNorm();
            const Vector2D& rf = bd.rFaceVec();

            divUStar(cell) += rho(bd.face())/timeStep*dot(u(bd.face()), sf) + (p(bd.face()) - p(cell))*dot(rf, sf)/dot(rf, rf);
        }
    }
}

void FractionalStep::computeAdvectingVelocity(Scalar timeStep)
{
    for(const Face &face: u.grid.interiorFaces())
    {
        const Cell &lCell = face.lCell();
        const Cell &rCell = face.rCell();

        Scalar g = rCell.volume()/(lCell.volume() + rCell.volume());
        Vector2D rc = rCell.centroid() - lCell.centroid();

        u(face) = g*(u(lCell) + timeStep*(gradP(lCell) - sg(lCell))/rho(lCell))
                + (1. - g)*(u(rCell) + timeStep*(gradP(rCell) - sg(rCell))/rho(rCell))
                - timeStep/rho(face)*(p(rCell) - p(lCell))*rc/rc.magSqr()
                + timeStep/rho(face)*g_;
    }

    for(const Face &face: u.grid.boundaryFaces())
    {
        const Cell &cell = face.lCell();

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT: case VectorFiniteVolumeField::OUTFLOW:
        {
            Vector2D rf = face.centroid() - cell.centroid();

            u(face) = u(cell)
                    + timeStep*(gradP(cell) - sg(cell))/rho(cell)
                    - timeStep/rho(face)*(p(face) - p(cell))*rf/rf.magSqr()
                    + timeStep/rho(face)*g_;
        }
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            Vector2D nWall = face.outwardNorm(face.lCell().centroid());
            u(face) = u(cell) - dot(u(cell), nWall)*nWall/nWall.magSqr();
            break;
        }

        default:
            throw Exception("FractionalStep", "computeAdvectingVelocity", "unrecongnized boundary condition type.");
        }
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

