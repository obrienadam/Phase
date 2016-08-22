#include "Piso.h"
#include "Exception.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"
#include "FaceInterpolation.h"
#include "GradientEvaluation.h"

Piso::Piso(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(addVectorField(input, "u")),
      h(addVectorField("h")),
      sg(addVectorField("sg")),
      gradP(addVectorField("gradP")),
      gradPCorr(addVectorField("gradPCorr")),
      p(addScalarField(input, "p")),
      pCorr(addScalarField("pCorr")),
      rho(addScalarField("rho")),
      mu(addScalarField("mu")),
      m(addScalarField("m")),
      d(addScalarField("d")),
      uEqn_(u, "momentum", SparseMatrix::IncompleteLUT),
      pCorrEqn_(pCorr, "pressure correction", SparseMatrix::IncompleteLUT)
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g", "(0,0)"));

    nInnerIterations_ = input.caseInput().get<size_t>("Solver.numInnerIterations");
    nPCorrections_ = input.caseInput().get<size_t>("Solver.numPressureCorrections");
    momentumOmega_ = input.caseInput().get<Scalar>("Solver.momentumRelaxation");
    pCorrOmega_ = input.caseInput().get<Scalar>("Solver.pressureCorrectionRelaxation");


    pCorr.copyBoundaryTypes(p);

    uEqn_.matrix().setFill(3);
    pCorrEqn_.matrix().setFill(3);

    //- Perform a pseudo initialization
    u.savePreviousTimeStep(0., 1);
}

Scalar Piso::solve(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    for(size_t innerIter = 0; innerIter < nInnerIterations_; ++innerIter)
    {
        u.savePreviousIteration();
        solveUEqn(timeStep);

        for(size_t pCorrIter = 0; pCorrIter < nPCorrections_; ++pCorrIter)
        {
            solvePCorrEqn();
            correctPressure();
            correctVelocity();
        }
    }

    printf("Max Co = %lf\n", courantNumber(timeStep));

    return 0.;
}

Scalar Piso::computeMaxTimeStep(Scalar maxCo) const
{
    Scalar maxTimeStepSqr = std::numeric_limits<Scalar>::infinity(), maxCoSqr = maxCo*maxCo;

    for(const Cell &cell: u.grid.fluidCells())
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaXSqr = nb.rCellVec().magSqr();
            Scalar magUSqr = u.faces()[nb.face().id()].magSqr();

            maxTimeStepSqr = std::min(maxTimeStepSqr, maxCoSqr*deltaXSqr/magUSqr);
        }

    return sqrt(maxTimeStepSqr);
}

//- Protected methods

Scalar Piso::solveUEqn(Scalar timeStep)
{
    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + cn::div(rho*u, u) + ib_.eqns(u) == ab::laplacian(mu, u) - fv::source(gradP) + fv::source(sg));
    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();

    return error;
}

Scalar Piso::solvePCorrEqn()
{
    for(const Cell& cell: m.grid.fluidCells())
    {
        Scalar &mass = m[cell.id()] = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            mass += dot(u.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            mass += dot(u.faces()[bd.face().id()], bd.outwardNorm());
    }

    pCorrEqn_ = (fv::laplacian(d, pCorr) + ib_.eqns(pCorr) == m);

    Scalar error = pCorrEqn_.solve();

    computeGradient(fv::GREEN_GAUSS_CELL_CENTERED, pCorr, gradPCorr);

    return error;
}

void Piso::rhieChowInterpolation()
{   
    VectorFiniteVolumeField& uStar = u.prevIter();
    VectorFiniteVolumeField& uPrev = u.prev(0);
    const Scalar dt = u.prevTimeStep(0);

    const auto diag = uEqn_.matrix().diagonal();

    d.fill(0.);
    for(const Cell& cell: d.grid.fluidCells())
        d[cell.id()] = cell.volume()/diag[cell.globalIndex()];

    interpolateFaces(fv::INVERSE_VOLUME, d);

    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell& cellP = face.lCell();
        const Cell& cellQ = face.rCell();

        const Vector2D rc = cellQ.centroid() - cellP.centroid();

        const Scalar dP = d(cellP);
        const Scalar dQ = d(cellQ);
        const Scalar df = d(face);
        const Scalar rhoP = rho(cellP);
        const Scalar rhoQ = rho(cellQ);
        const Scalar rhof = rho(face);

        const Scalar g = cellQ.volume()/(cellP.volume() + cellQ.volume());

        if(ib_.isIbCell(cellP) || ib_.isIbCell(cellQ))
        {
            u(face) = g*u(cellP) + (1. - g)*u(cellQ);
        }
        else
        {
            u(face) = g*u(cellP) + (1. - g)*u(cellQ)
                    + (1. - momentumOmega_)*(uStar(face) - (g*uStar(cellP) + (1. - g)*uStar(cellQ)))
                    + (rhof*df*uPrev(face) - (g*rhoP*dP*uPrev(cellP) + (1. - g)*rhoQ*dQ*uPrev(cellQ)))/dt //- This term is very important!
                    - df*(p(cellQ) - p(cellP))*rc/dot(rc, rc) + rhof*(g*dP*gradP(cellP)/rhoP + (1. - g)*dQ*gradP(cellQ)/rhoQ)/2.
                    + df*rhof*g_ - rhof*(g*dP*sg(cellP)/rhoP + (1. - g)*dQ*sg(cellQ)/rhoQ)/2.;
        }
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Vector2D rf = face.centroid() - face.lCell().centroid();
        const Scalar df = d(face);
        const Scalar rhof = rho(face);

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) = u(face.lCell())
                    - df*((p(face) - p(face.lCell()))*rf/dot(rf, rf) - rhof*gradP(face.lCell())/rho(face.lCell()))
                    + df*(rhof*g_ - rhof*sg(face.lCell())/rho(face.lCell()));
            break;

        case VectorFiniteVolumeField::SYMMETRY:
        {
            const Vector2D nWall = face.outwardNorm(face.lCell().centroid());

            u(face) = u(face.lCell()) - dot(u(face.lCell()), nWall)*nWall/nWall.magSqr();
        }
            break;

        case VectorFiniteVolumeField::OUTFLOW:
            u(face) = u(face.lCell())
                    - df*((p(face) - p(face.lCell()))*rf/dot(rf, rf) - rhof*gradP(face.lCell())/rho(face.lCell()))
                    + df*(rhof*g_ - rhof*sg(face.lCell())/rho(face.lCell()));

            if(dot(u(face), face.outwardNorm(face.lCell().centroid())) < 0.)
                u(face) = Vector2D(0., 0.);
            break;

        default:
            throw Exception("Piso", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}

void Piso::correctPressure()
{
    for(const Cell& cell: p.grid.activeCells())
        p(cell) += pCorrOmega_*pCorr(cell);

    interpolateFaces(fv::INVERSE_VOLUME, p);
    computeStaticPressure();
    computeGradient(fv::GREEN_GAUSS_CELL_CENTERED, p, gradP);
}

void Piso::correctVelocity()
{
    for(const Cell& cell: u.grid.fluidCells())
        u(cell) -= d(cell)*gradPCorr(cell);

    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        Vector2D rc = rCell.centroid() - lCell.centroid();

        u(face) -= d(face)*(pCorr(rCell) - pCorr(lCell))*rc/dot(rc, rc);
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Vector2D rf = face.centroid() - face.lCell().centroid();

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;
        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u(face) -= d(face)*(pCorr(face) - pCorr(face.lCell()))*rf/dot(rf, rf);
            break;
        }
    }

    //- Update mass source for the purpose of error checking

    for(const Cell& cell: m.grid.fluidCells())
    {
        m(cell) = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            m(cell) += dot(u(nb.face()), nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            m(cell) += dot(u(bd.face()), bd.outwardNorm());
    }
}

void Piso::computeStaticPressure()
{
    for(const Face& face: grid_.boundaryFaces())
    {
        if(p.boundaryType(face.id()) == ScalarFiniteVolumeField::NORMAL_GRADIENT)
            p(face) = p(face.lCell()) + rho(face.lCell())*dot(g_, face.centroid() - face.lCell().centroid());
    }

    for(const ImmersedBoundaryObject& ibObj: ib_.ibObjs())
        for(const Cell &cell: ibObj.cells())
            for(const InteriorLink &nb: cell.neighbours())
            {
                if(grid_.fluidCells().isInGroup(nb.cell()) && ibObj.boundaryType(p.name) == ImmersedBoundaryObject::NORMAL_GRADIENT)
                {
                    p(nb.face()) = p(nb.cell()) + rho(nb.cell())*dot(g_, nb.face().centroid() - nb.cell().centroid());
                }
            }
}

Scalar Piso::courantNumber(Scalar timeStep)
{
    Scalar maxCoSqr = 0., timeStepSqr = timeStep*timeStep;

    for(const Cell &cell: u.grid.fluidCells())
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaXSqr = nb.rCellVec().magSqr();
            Scalar magUSqr = u.faces()[nb.face().id()].magSqr();

            maxCoSqr = std::max(maxCoSqr, magUSqr*timeStepSqr/deltaXSqr);
        }

    return sqrt(maxCoSqr);
}
