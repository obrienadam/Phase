#include "Plic.h"

namespace plic
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field, Scalar timeStep, std::vector<Polygon> &plicPolygons)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);
    VectorFiniteVolumeField gradField = grad(field);

    entries.reserve(field.grid.nActiveCells());

    for(const Cell &cell: field.grid.fluidCells())
    {
        Scalar centralCoeff = 0.;
        const size_t row = cell.globalIndex();

        Polygon plicPgn = computeInterfacePolygon(cell, field[cell.id()], -gradField[cell.id()]);

        //plicPolygons[cell.id()] = plicPgn;

        for(const InteriorLink &nb: cell.neighbours())
        {
            const size_t col = nb.cell().globalIndex();

            if(dot(u.faces()[nb.face().id()], nb.outwardNorm()) > 0.)
            {
                Polygon fluxPgn = computeFluxPolygon(nb, u.faces()[nb.face().id()], timeStep);

                if(fluxPgn.area() == 0.)
                    continue;

                fluxPgn = intersectionPolygon(plicPgn, fluxPgn);

                //if(fluxPgn.area() > 1.01*plicPgn.area())
                //   throw Exception("plic", "div", "invalid flux polygon. Face flux = " + std::to_string(fluxPgn.area()) + ", flux polygon area = " + std::to_string(plicPgn.area()));

                eqn.boundaries()(row) -= fluxPgn.area();
                eqn.boundaries()(col) += fluxPgn.area();
            }
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {
            Polygon fluxPgn = computeFluxPolygon(bd, u.faces()[bd.face().id()], timeStep);

            switch(field.boundaryType(bd.face().id()))
            {
            case ScalarFiniteVolumeField::FIXED:
                eqn.boundaries()(row) -= fluxPgn.area()*field.faces()[bd.face().id()]; // No intersection computation required, always same gamma influx
                break;

            case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                fluxPgn = intersectionPolygon(plicPgn, fluxPgn); // Flux depends on the plic reconstruction in the boundary cell
                eqn.boundaries()(row) -= fluxPgn.area();
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
    if(gamma >= 1 - 1e-10)
        return cell.shape();
    else if (gamma <= 1e-10)
        return Polygon();
    else if (interfaceNormal.x == 0 && interfaceNormal.y == 0) // This is bad and should not happen, interface cells should have a gradient. Plic pgn will be placed in the middle4
    {
        Scalar a = 0., b = 1.;

        while(true)
        {
            Scalar c = (a + b)/2.;

            Polygon plicPgn = cell.shape();
            plicPgn.scale(c);

            Scalar f = plicPgn.area()/cell.volume() - gamma;

            if(fabs(f) < 1e-12)
                return plicPgn;
            else if(f > 0)
                b = c;
            else
                a = c;
        }
    }

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

    if(fabs(plicPgn.area()/cell.volume() - gamma) > 1e-5)
        throw Exception("plic", "computeInterfacePolygon", "an invalid PLIC polygon reconstruction was detected. Reconstructed gamma = "
                        + std::to_string(plicPgn.area()/cell.volume()) + ", actual gamma = " + std::to_string(gamma) + ".");

    return plicPgn;
}

Polygon computeFluxPolygon(const BoundaryLink &link, const Vector2D& uf, Scalar timeStep)
{
    Vector2D off = dot(uf, link.outwardNorm())*link.outwardNorm()/link.outwardNorm().magSqr()*timeStep;

    std::vector<Point2D> verts = {
        link.face().lNode(),
        link.face().lNode() - off,
        link.face().rNode() - off,
        link.face().rNode(),
    };

    return Polygon(verts);
}
}
