#include "Plic.h"

namespace plic
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field, Scalar timeStep)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);
    VectorFiniteVolumeField gradField = grad(field); //- should this be smoothed?

    entries.reserve(field.grid.nActiveCells());

    // Figure out the plic reconstruction
    for(const Cell &cell: field.grid.fluidCells())
    {
        bool isInterfaceCell = gradField[cell.id()].magSqr() > 1e-8;

        Scalar centralCoeff = 0.;
        const size_t row = cell.globalIndex();

        for(const InteriorLink &nb: cell.neighbours())
        {
            const size_t col = nb.cell().globalIndex();

            if(dot(u.faces()[nb.face().id()], nb.outwardNorm()) >= 0.)
            {
                //- Get the plic polygon
                Polygon plicPgn = computeInterfacePolygon(cell, field[cell.id()], -gradField[cell.id()]);

                //- Get the flux polygon
                Polygon fluxPgn = computeFluxPolygon(nb, u.faces()[nb.face().id()], timeStep);

                //- Compute the flux
                Scalar faceFlux = intersectionPolygon(plicPgn, fluxPgn).area();

                eqn.boundaries()(row) -= faceFlux;
                eqn.boundaries()(col) += faceFlux;
            }
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Polygon fluxPgn = computeFluxPolygon(bd, u.faces()[bd.face().id()], timeStep);
            Polygon plicPgn;
            Vector2D wt, n;
            Scalar faceFlux;

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= fluxPgn.area()*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                wt = bd.outwardNorm().tangentVec();
                n = dot(wt, gradField[cell.id()]) > 0 ? wt : -wt;

                plicPgn = computeInterfacePolygon(cell, field[cell.id()], n);
                fluxPgn = computeFluxPolygon(bd, u.faces()[bd.face().id()], timeStep);
                faceFlux = intersectionPolygon(plicPgn, fluxPgn).area();

                eqn.boundaries()(row) -= faceFlux;

                break;

            default:
                throw Exception("plic", "div", "unrecognized or unspecified boundary type.");
            }
        }

        centralCoeff += 1.;
        eqn.boundaries()(row) += field[cell.id()];

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Polygon computeInterfacePolygon(const Cell &cell, Scalar gamma,  const Vector2D& interfaceNormal)
{
    if(gamma >= 1)
        return cell.shape();
    else if (gamma <= 0)
        return Polygon();

    Line2D plic(cell.centroid(), interfaceNormal);

    auto f = [&cell, &plic, &gamma](Scalar c)
    {
        return clipPolygon(cell.shape(), plic.adjust(c)).area()/cell.volume() - gamma;
    };

    auto bracket = [&cell, &plic, &gamma](Scalar &cMin, Scalar &cMax)
    {
        bool init = true;

        for(const Point2D& vtx: cell.shape())
        {
            Scalar c = dot(vtx - plic.r0(), plic.n());

            if(init)
            {
                cMin = cMax = c;
                init = false;
            }
            else
            {
                cMin = c < cMin ? c : cMin;
                cMax = c > cMax ? c : cMax;
            }

        }
    };

    Scalar func, c, cMin, cMax;
    const int maxIters = 50;
    int iters = 0;
    //- find the bounds
    bracket(cMin, cMax);
    //- Make a guess

    //- bisection
    while(true)
    {
        c = (cMin + cMax)/2;
        func = f(c);

        if (fabs(func) < 1e-12)
            break;
        else if (++iters > maxIters)
            break;
        else if (func < 0)
            cMin = c;
        else
            cMax = c;
    }

    return clipPolygon(cell.shape(), plic.adjust(c)); // Return the clipped polygon that represents the shape of phase B
}

Polygon computeFluxPolygon(const BoundaryLink &link, const Vector2D& uf, Scalar timeStep)
{
    Vector2D off = -uf*timeStep;

    std::vector<Point2D> verts = {
        link.face().lNode(),
        link.face().lNode() + off,
        link.face().rNode() + off,
        link.face().rNode(),
    };

    return Polygon(verts);
}
}
