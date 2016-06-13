#include "Simple.h"
#include "Exception.h"

Simple::Simple(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(addVectorField(input, "u")),
      h(addVectorField("h")),
      sg(addVectorField("sg")),
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
    momentumOmega_ = input.caseInput().get<Scalar>("Solver.momentumRelaxation");
    pCorrOmega_ = input.caseInput().get<Scalar>("Solver.pressureCorrectionRelaxation");


    pCorr.copyBoundaryTypes(p);
}

Scalar Simple::solve(Scalar timeStep)
{
    u.savePreviousTimeStep(timeStep, 1);

    for(size_t i = 0; i < nInnerIterations_; ++i)
    {
        u.savePreviousIteration();

        solveUEqn(timeStep);
        solvePCorrEqn();
        correctPressure();
        correctVelocity();
    }

    printf("Max Co = %lf\n", courantNumber(timeStep));

    return 0.;
}

Scalar Simple::computeMaxTimeStep(Scalar maxCo) const
{
    Scalar maxTimeStep = std::numeric_limits<Scalar>::infinity();

    for(const Cell &cell: u.grid.fluidCells())
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaX = nb.rCellVec().mag();
            Scalar magU = u.faces()[nb.face().id()].mag();

            maxTimeStep = std::min(maxTimeStep, maxCo*deltaX/magU);
        }

    return maxTimeStep;
}

//- Protected methods

Scalar Simple::solveUEqn(Scalar timeStep)
{
    sg = fv::gravity(rho, g_);

    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) == fv::laplacian(mu, u) - fv::grad(p) + fv::source(sg));
    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();

    return error;
}

Scalar Simple::solvePCorrEqn()
{
    for(const Cell& cell: m.grid.fluidCells())
    {
        size_t id = cell.id();

        m[id] = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            m[id] += dot(u.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            m[id] += dot(u.faces()[bd.face().id()], bd.outwardNorm());
    }

    pCorrEqn_ = (fv::laplacian(d, pCorr) == m);

    Scalar error = pCorrEqn_.solve();

    interpolateFaces(pCorr);

    return error;
}

void Simple::computeD()
{
    const auto diag = uEqn_.matrix().diagonal();

    d.fill(0.);
    for(const Cell& cell: d.grid.fluidCells())
        d[cell.id()] = cell.volume()/diag[cell.globalIndex()];

    interpolateFaces(d);
}

void Simple::rhieChowInterpolation()
{
    computeD();
    VectorFiniteVolumeField gradP = grad(p);
    const VectorFiniteVolumeField& uStar = u.prevIter();

    for(const Face& face: u.grid.interiorFaces())
    {
        const Cell& cellP = face.lCell();
        const Cell& cellQ = face.rCell();
        const size_t fid = face.id();

        const Scalar g = cellQ.volume()/(cellP.volume() + cellQ.volume());
        const Vector2D rc = cellQ.centroid() - cellP.centroid();

        const Scalar df = d.faces()[fid];
        const Scalar rhoP = rho[cellP.id()];
        const Scalar rhoQ = rho[cellQ.id()];
        const Scalar rhof = rho.faces()[fid];

        u.faces()[fid] = g*u[cellP.id()] + (1. - g)*u[cellQ.id()]
                + (1. - momentumOmega_)*(uStar.faces()[fid] - (g*uStar[cellP.id()] + (1. - g)*uStar[cellQ.id()]))
                - df*((p[cellQ.id()] - p[cellP.id()])*rc/dot(rc, rc) - rhof*(gradP[cellP.id()]/rhoP + gradP[cellQ.id()]/rhoQ)/2.)
                + df*(rhof*g_ - rhof*(sg[cellP.id()]/rhoP + sg[cellQ.id()]/rhoQ)/2.);
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
                    - df*((p.faces()[fid] - p[cellP.id()])*rf/dot(rf, rf) - gradP[cellP.id()])
                    + df*(rhof*g_ - sg[cellP.id()]);
            break;

        case VectorFiniteVolumeField::SYMMETRY:
            nWall = face.outwardNorm(face.lCell().centroid()).unitVec();
            u.faces()[face.id()] = u[face.lCell().id()] - dot(u[face.lCell().id()], nWall)*nWall;
            break;

        default:
            throw Exception("Simple", "rhieChowInterpolation", "unrecognized boundary condition type.");
        }
    }
}

void Simple::correctPressure()
{
    for(const Cell& cell: p.grid.activeCells())
        p[cell.id()] += pCorrOmega_*pCorr[cell.id()];

    interpolateFaces(p);
}

void Simple::correctVelocity()
{
    VectorFiniteVolumeField gradPCorr = grad(pCorr);

    for(const Cell& cell: u.grid.fluidCells())
        u[cell.id()] -= d[cell.id()]*gradPCorr[cell.id()];

    for(const Face& face: u.grid.interiorFaces())
    {
        const size_t faceId = face.id();

        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        Vector2D sf = face.outwardNorm(lCell.centroid());
        Vector2D rc = rCell.centroid() - lCell.centroid();

        u.faces()[faceId] -= d.faces()[faceId]*(pCorr[rCell.id()] - pCorr[lCell.id()])*sf/dot(rc, rc);
    }

    for(const Face& face: u.grid.boundaryFaces())
    {
        switch(u.boundaryType(face.id()))
        {
        case VectorFiniteVolumeField::FIXED:
            break;
        case VectorFiniteVolumeField::NORMAL_GRADIENT:
            u.faces()[face.id()] = u[face.lCell().id()];
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

Scalar Simple::courantNumber(Scalar timeStep)
{
    Scalar maxCo = 0.;

    for(const Cell &cell: u.grid.fluidCells())
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaX = nb.rCellVec().mag();
            Scalar magU = u.faces()[nb.face().id()].mag();

            maxCo = std::max(maxCo, magU*timeStep/deltaX);
        }

    return maxCo;
}
