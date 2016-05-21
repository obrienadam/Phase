#include "Plic.h"

namespace plic
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field, Scalar timeStep, std::vector<Polygon> &plicPolygons)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);
    VectorFiniteVolumeField gradField = grad(field);

    entries.reserve(field.grid.nActiveCells());

    // Figure out the plic reconstruction
    for(const Cell &cell: field.grid.fluidCells())
    {
        Scalar centralCoeff = 0.;
        const size_t row = cell.globalIndex();

        //- Get the plic polygon
        Polygon &plicPgn = plicPolygons[cell.id()] = computeInterfacePolygon(cell, field[cell.id()], -gradField[cell.id()]);

        for(const InteriorLink &nb: cell.neighbours())
        {
            const size_t col = nb.cell().globalIndex();

            if(dot(u.faces()[nb.face().id()], nb.outwardNorm()) > 0.)
            {
                //- Get the flux polygon
                Polygon fluxPgn = computeFluxPolygon(nb, u.faces()[nb.face().id()], timeStep);

                //- Compute the flux
                Scalar faceFlux = plicPgn.area()*0.001;

                eqn.boundaries()(row) -= faceFlux;
                //eqn.boundaries()(col) += faceFlux;
            }
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Polygon fluxPgn = computeFluxPolygon(bd, u.faces()[bd.face().id()], timeStep);
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

                eqn.boundaries()(row) -= plicPgn.area()*0.001;

                break;

            default:
                throw Exception("plic", "div", "unrecognized or unspecified boundary type.");
            }
        }

        centralCoeff += cell.volume();
        eqn.boundaries()(row) += field[cell.id()]*cell.volume();

        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff));
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Polygon computeInterfacePolygon(const Cell &cell, Scalar gamma,  const Vector2D& interfaceNormal)
{
    if(gamma >= 1 - 1e-8)
        return cell.shape();
    else if (gamma <= 1e-8)
        return Polygon();
    else if (interfaceNormal.x == 0 && interfaceNormal.y == 0)
        return cell.shape();

    Line2D plic(cell.centroid(), interfaceNormal);

    auto f = [&cell, &plic, &gamma](Scalar c)->Scalar
    {
        return clipPolygon(cell.shape(), plic.adjust(c)).area()/cell.volume() - gamma;
    };

    auto bracket = [&cell, &plic, &gamma](Scalar &cMin, Scalar &cMax)->void
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

    Polygon plicPgn = clipPolygon(cell.shape(), plic.adjust(c)); // Return the clipped polygon that represents the shape of phase B

    if(fabs(plicPgn.area()/cell.volume() - gamma) > 1e-2)
        throw Exception("plic", "computeInterfacePolygon", "an invalid PLIC polygon reconstruction was detected. Reconstructed gamma = "
                        + std::to_string(plicPgn.area()/cell.volume()) + ", actual gamma = " + std::to_string(gamma) + ".");

    return plicPgn;
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
