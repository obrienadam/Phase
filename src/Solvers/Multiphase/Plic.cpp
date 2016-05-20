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

            if(isInterfaceCell)
            {
                if(dot(u.faces()[nb.face().id()], nb.outwardNorm()) > 0.)
                {
                    //- Get the plic polygon
                    Polygon plicPgn = computeInterfacePolygon(cell, field[cell.id()], -gradField[cell.id()]);

                    //- Get the flux polygon
                    Polygon fluxPgn = computeFluxPolygon(nb, u.faces()[nb.face().id()], timeStep);

                    //- Compute the flux
                    Scalar faceFlux = intersectionPolygon(plicPgn, fluxPgn).area()/(timeStep);

                    eqn.boundaries()(row) -= faceFlux;
                    eqn.boundaries()(col) += faceFlux;
                }
            }
            else
            {
                Scalar faceFlux = std::max(dot(u.faces()[nb.face().id()], nb.outwardNorm()), 0.);

                eqn.boundaries()(row) -= faceFlux*field[cell.id()];
                eqn.boundaries()(col) += faceFlux*field[cell.id()];
            }
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Scalar faceFlux = dot(u.faces()[bd.face().id()], bd.outwardNorm());

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= faceFlux*field.faces()[bd.face().id()];
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                centralCoeff += faceFlux;
                eqn.boundaries()(row) -= bd.outwardNorm().mag()*field.boundaryRefValue(bd.face().id());
                break;

            default:
                throw Exception("plic", "div", "unrecognized or unspecified boundary type.");
            }
        }

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Polygon computeInterfacePolygon(const Cell &cell, Scalar gamma,  const Vector2D& interfaceNormal)
{
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

Polygon computeFluxPolygon(const InteriorLink &link, const Vector2D& uf, Scalar timeStep)
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
