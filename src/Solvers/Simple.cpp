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
      uEqn_(grid, u),
      pCorrEqn_(grid, pCorr)
{
    rho.fill(input.caseInput().get<Scalar>("Properties.rho"));
    mu.fill(input.caseInput().get<Scalar>("Properties.mu"));
    g_ = Vector2D(input.caseInput().get<std::string>("Properties.g"));
}

Scalar Simple::solve(Scalar timeStep)
{
    solveUEqn(timeStep);
    //solvePCorrEqn();
    return 0.;
}

//- Protected methods

Scalar Simple::solveUEqn(Scalar timeStep)
{
    interpolateFaces(p);

    uEqn_ = (rho*ddt(u, timeStep) + rho*div(u, u) == mu*laplacian(u) - grad(p));
    return uEqn_.solve();
}

Scalar Simple::solvePCorrEqn()
{
    rhieChowInterpolation();

    pCorrEqn_ = (rho*d*laplacian(pCorr) == m);
    return pCorrEqn_.solve();
}

void Simple::computeD()
{
    const auto diag = uEqn_.matrix().diagonal();

    for(const Cell& cell: d.grid.cells)
        d[cell.id()] = cell.volume()/diag[cell.id()];
}

void Simple::rhieChowInterpolation()
{
    computeD();

    h = u + d*grad(p);

    interpolateFaces(h);
    interpolateFaces(d);

    for(const Face& face: u.grid.faces)
    {
        if(face.isBoundary())
            continue;

        const Cell& lCell = face.lCell();
        const Cell& rCell = face.rCell();

        Vector2D sf = face.outwardNorm(lCell.centroid());
        Vector2D rc = rCell.centroid() - lCell.centroid();

        size_t faceId = face.id();

        u.faces()[faceId] = h.faces()[faceId] - d.faces()[faceId]*(p[rCell.id()] - p[lCell.id()])*sf/dot(rc, rc);
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
