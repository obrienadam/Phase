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
      uEqn_(u),
      pCorrEqn_(pCorr)
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho"));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu"));
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g"));

    momentumOmega_ = input.caseInput().get<Scalar>("Solver.momentumRelaxation");
    pCorrOmega_ = input.caseInput().get<Scalar>("Solver.pressureCorrectionRelaxation");

    pCorr.copyBoundaryTypes(p);
}

Scalar Simple::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    solvePCorrEqn();
    correctPressure();
    //correctVelocity();
    return 0.;
}

//- Protected methods

Scalar Simple::solveUEqn(Scalar timeStep)
{
    interpolateFaces(p);

    uEqn_ = (fv::div(rho*u, u) == fv::laplacian(mu, u) - fv::grad(p));
    uEqn_.relax(momentumOmega_);

    return uEqn_.solve();
}

Scalar Simple::solvePCorrEqn()
{
    rhieChowInterpolation();

    pCorrEqn_ = (fv::laplacian(rho*d, pCorr) == m);
    return pCorrEqn_.solve();
}

void Simple::computeD()
{
    const auto diag = uEqn_.matrix().diagonal();

    for(const Cell& cell: d.grid.cells)
        d[cell.id()] = cell.volume()/diag[cell.id()];

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
        else
        {
            u.faces()[faceId] = h.faces()[faceId];
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
}

void Simple::correctVelocity()
{
    interpolateFaces(pCorr);
    VectorFiniteVolumeField gradPCorr = grad(pCorr);

    for(const Cell& cell: u.grid.cells)
    {
        if(!cell.isActive())
            continue;

        u[cell.id()] -= d[cell.id()]*gradPCorr[cell.id()];
    }

    for(const Face& face: u.grid.faces)
    {
        if(face.isBoundary())
            continue;

        const size_t faceId = face.id();

        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        Vector2D sf = face.outwardNorm(lCell.centroid());
        Vector2D rc = rCell.centroid() - lCell.centroid();

        u.faces()[faceId] += d.faces()[faceId]*(pCorr[rCell.id()] - pCorr[lCell.id()])*sf/dot(rc, rc);
    }

    //rhieChowInterpolation();
}
