#include "Piso.h"
#include "Exception.h"
#include "CrankNicolson.h"
#include "AdamsBashforth.h"
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

    interpolateFaces(d);

    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell& cellP = face.lCell();
        const Cell& cellQ = face.rCell();
        const size_t fid = face.id();

        const Vector2D rc = cellQ.centroid() - cellP.centroid();

        const Scalar dP = d[cellP.id()];
        const Scalar dQ = d[cellQ.id()];
        const Scalar df = d.faces()[fid];
        const Scalar rhoP = rho[cellP.id()];
        const Scalar rhoQ = rho[cellQ.id()];
        const Scalar rhof = rho.faces()[fid];

        const Scalar g = cellQ.volume()/(cellP.volume() + cellQ.volume());

        u.faces()[fid] = g*u[cellP.id()] + (1. - g)*u[cellQ.id()]
                + (1. - momentumOmega_)*(uStar.faces()[fid] - (g*uStar[cellP.id()] + (1. - g)*uStar[cellQ.id()]))
                + (rhof*df*uPrev.faces()[fid] - (g*rhoP*dP*uPrev[cellP.id()] + (1. - g)*rhoQ*dQ*uPrev[cellQ.id()]))/dt //- This term is very important!
                - df*(p[cellQ.id()] - p[cellP.id()])*rc/dot(rc, rc) + rhof*(g*dP*gradP[cellP.id()]/rhoP + (1. - g)*dQ*gradP[cellQ.id()]/rhoQ)/2.
                + df*rhof*g_ - rhof*(g*dP*sg[cellP.id()]/rhoP + (1. - g)*dQ*sg[cellQ.id()]/rhoQ)/2.;
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        Vector2D nWall;

        const Cell& cellP = face.lCell();
        const size_t fid = face.id();
        const Vector2D rf = face.centroid() - cellP.centroid();
        const Scalar df = d.faces()[fid];
        const Scalar rhof = rho.faces()[fid];

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;

        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u.faces()[face.id()] = u[face.lCell().id()]
                    - df*((p.faces()[fid] - p[cellP.id()])*rf/dot(rf, rf) - rhof*gradP[cellP.id()]/rho[cellP.id()])
                    + df*(rhof*g_ - rhof*sg[cellP.id()]/rho[cellP.id()]);
            break;

        case VectorFiniteVolumeField::SYMMETRY:
            nWall = face.outwardNorm(face.lCell().centroid()).unitVec();
            u.faces()[face.id()] = u[face.lCell().id()] - dot(u[face.lCell().id()], nWall)*nWall;
            break;

        case VectorFiniteVolumeField::OUTFLOW:
            u.faces()[face.id()] = u[face.lCell().id()]
                    - df*((p.faces()[fid] - p[cellP.id()])*rf/dot(rf, rf) - rhof*gradP[cellP.id()]/rho[cellP.id()])
                    + df*(rhof*g_ - rhof*sg[cellP.id()]/rho[cellP.id()]);

            if(dot(u.faces()[face.id()], face.outwardNorm(cellP.centroid())) < 0.)
                u.faces()[face.id()] = Vector2D(0., 0.);
            break;

        default:
            throw Exception("Piso", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}

void Piso::correctPressure()
{
    for(const Cell& cell: p.grid.activeCells())
        p[cell.id()] += pCorrOmega_*pCorr[cell.id()];

    computeGradient(fv::GREEN_GAUSS_CELL_CENTERED, p, gradP);
    computeStaticPressure();

    //for(const Face &face: p.grid.boundaryFaces())
    //    p.faces()[face.id()] += rho[face.lCell().id()]*dot(face.centroid() - face.lCell().centroid(), g_);
}

void Piso::correctVelocity()
{
    for(const Cell& cell: u.grid.fluidCells())
        u[cell.id()] -= d[cell.id()]*gradPCorr[cell.id()];

    for(const Face& face: u.grid.interiorFaces())
    {
        const size_t faceId = face.id();

        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        Vector2D rc = rCell.centroid() - lCell.centroid();

        u.faces()[faceId] -= d.faces()[faceId]*(pCorr[rCell.id()] - pCorr[lCell.id()])*rc/dot(rc, rc);
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        const Vector2D rf = face.centroid() - face.lCell().centroid();

        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;
        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u.faces()[face.id()] -= d.faces()[face.id()]*(pCorr.faces()[face.id()] - pCorr[face.lCell().id()])*rf/dot(rf, rf);
            break;
        }
    }

    //- Update mass source for the purpose of error checking

    for(const Cell& cell: m.grid.fluidCells())
    {
        Scalar &mDot = m[cell.id()] = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            mDot += dot(u.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            mDot += dot(u.faces()[bd.face().id()], bd.outwardNorm());
    }
}

void Piso::computeStaticPressure()
{
    for(const Face& face: grid_.boundaryFaces())
    {
        if(p.boundaryType(face.id()) == ScalarFiniteVolumeField::NORMAL_GRADIENT)
            p.faces()[face.id()] = p[face.lCell().id()] + rho[face.lCell().id()]*dot(g_, face.centroid() - face.lCell().centroid());
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
