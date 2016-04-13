#include "Simple.h"
#include "Exception.h"

Simple::Simple(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      Solver(grid, input),
      u(grid.addVectorField(input, "u")),
      h(grid.addVectorField("h")),
      p(grid.addScalarField(input, "p")),
      pCorr(grid.addScalarField("pCorr")),
      rho(grid.addScalarField("rho")),
      mu(grid.addScalarField("mu")),
      m(grid.addScalarField("m")),
      d(grid.addScalarField("d")),
      ids(grid.addScalarField("ids")),
      globalIndices(grid.addScalarField("globalIndices")),
      uEqn_(u, "momentum"),
      pCorrEqn_(pCorr, "pressure correction")
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho", 1.));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu", 1.));
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g", "(0,0)"));

    nInnerIterations_ = input.caseInput().get<size_t>("Solver.numInnerIterations");
    momentumOmega_ = input.caseInput().get<Scalar>("Solver.momentumRelaxation");
    pCorrOmega_ = input.caseInput().get<Scalar>("Solver.pressureCorrectionRelaxation");


    pCorr.copyBoundaryTypes(p);

    for(const Cell& cell: ids.grid.cells)
    {
        ids[cell.id()] = cell.id();
        globalIndices[cell.id()] = cell.globalIndex();
    }
}

Scalar Simple::solve(Scalar timeStep)
{
    u.save();

    for(size_t i = 0; i < nInnerIterations_; ++i)
    {
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

    for(const Cell &cell: u.grid.cells)
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
    uEqn_ = (fv::ddt(rho, u, timeStep) + fv::div(rho*u, u) == fv::laplacian(mu, u) - fv::grad(p));
    uEqn_.relax(momentumOmega_);

    Scalar error = uEqn_.solve();

    rhieChowInterpolation();

    return error;
}

Scalar Simple::solvePCorrEqn()
{
    pCorrEqn_ = (fv::laplacian(rho*d, pCorr) == m);

    Scalar error = pCorrEqn_.solve();

    interpolateFaces(pCorr);

    return error;
}

void Simple::computeD()
{
    const auto diag = uEqn_.matrix().diagonal();

    for(const Cell& cell: d.grid.cells)
    {
        if(cell.isActive())
            d[cell.id()] = cell.volume()/diag[cell.globalIndex()];
    }

    interpolateFaces(d);
}

void Simple::rhieChowInterpolation()
{
    computeD();

    h = u + d*grad(p);
    interpolateFaces(h);

    for(const Face& face: u.grid.faces)
    {
        size_t faceId = face.id();

        if(!face.isBoundary())
        {
            const Cell& lCell = face.lCell();
            const Cell& rCell = face.rCell();

            Vector2D sf = face.outwardNorm(lCell.centroid());
            Vector2D rc = rCell.centroid() - lCell.centroid();

            u.faces()[faceId] = h.faces()[faceId] - d.faces()[faceId]*(p[rCell.id()] - p[lCell.id()])*sf/dot(rc, rc);
        }
        else if(u.boundaryType(faceId) == VectorFiniteVolumeField::NORMAL_GRADIENT)
        {
            u.faces()[faceId] = u[face.lCell().id()];
        }
    }

    for(const Cell& cell: m.grid.cells)
    {
        size_t id = cell.id();

        m[id] = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            m[id] += rho.faces()[nb.face().id()]*dot(u.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            m[id] += rho.faces()[bd.face().id()]*dot(u.faces()[bd.face().id()], bd.outwardNorm());
    }
}

void Simple::correctPressure()
{
    for(const Cell& cell: p.grid.cells)
        p[cell.id()] += pCorrOmega_*pCorr[cell.id()];

    interpolateFaces(p);
}

void Simple::correctVelocity()
{
    VectorFiniteVolumeField gradPCorr = grad(pCorr);

    for(const Cell& cell: u.grid.cells)
    {
        if(!cell.isActive())
            continue;

        u[cell.id()] -= d[cell.id()]*gradPCorr[cell.id()];
    }
    /*
    for(const Face& face: u.grid.faces)
    {
        const size_t faceId = face.id();

        if(face.isInterior())
        {
            const Cell& lCell = face.lCell();
            const Cell& rCell = face.rCell();

            Vector2D sf = face.outwardNorm(lCell.centroid());
            Vector2D rc = rCell.centroid() - lCell.centroid();

            u.faces()[faceId] += d.faces()[faceId]*(pCorr[rCell.id()] - pCorr[lCell.id()])*sf/dot(rc, rc);
        }
        else if(u.boundaryType(faceId) == VectorFiniteVolumeField::NORMAL_GRADIENT)
        {
            u.faces()[faceId] = u[face.lCell().id()];
        }
    }
    */

    rhieChowInterpolation(); // This can be used instead of the velocity face correction, but is probably more expensive

    //- Update mass source for the purpose of error checking

    for(const Cell& cell: m.grid.cells)
    {
        Scalar &mDot = m[cell.id()] = 0.;

        for(const InteriorLink& nb: cell.neighbours())
            mDot += rho.faces()[nb.face().id()]*dot(u.faces()[nb.face().id()], nb.outwardNorm());

        for(const BoundaryLink& bd: cell.boundaries())
            mDot += rho.faces()[bd.face().id()]*dot(u.faces()[bd.face().id()], bd.outwardNorm());
    }
}

Scalar Simple::courantNumber(Scalar timeStep)
{
    Scalar maxCo = 0.;

    for(const Cell &cell: u.grid.cells)
        for(const InteriorLink &nb: cell.neighbours())
        {
            Scalar deltaX = nb.rCellVec().mag();
            Scalar magU = u.faces()[nb.face().id()].mag();

            maxCo = std::max(maxCo, magU*timeStep/deltaX);
        }

    return maxCo;
}
