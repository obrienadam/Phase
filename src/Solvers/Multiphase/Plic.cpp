#include "Plic.h"

namespace plic
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field, Scalar timeStep, std::vector<Polygon> &plicPolygons)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;

    for(int componentNo = 0; componentNo < 2; ++componentNo)
    {
        Equation<ScalarFiniteVolumeField> eqn(field);
        VectorFiniteVolumeField gradField = grad(field);
        entries.clear();
        entries.reserve(field.grid.nActiveCells());

        for(const Cell &cell: field.grid.fluidCells())
        {
            const size_t row = cell.globalIndex();

            Polygon plicPgn = computeInterfacePolygon(cell, field[cell.id()], -gradField[cell.id()]);

            for(const InteriorLink &nb: cell.neighbours())
            {
                const size_t col = nb.cell().globalIndex();

                if(dot(u.faces()[nb.face().id()], nb.outwardNorm()) > 0.)
                {
                    Polygon fluxPgn = computeFluxPolygon(nb, u.faces()[nb.face().id()], timeStep, componentNo);

                    if(fluxPgn.isEmpty())
                        continue;

                    Scalar massTransfer;

                    if(!plicPgn.isEmpty())
                        massTransfer = intersectionPolygon(plicPgn, fluxPgn).area();
                    else
                        massTransfer = fluxPgn.area()*field[cell.id()];

                    eqn.boundaries()(row) -= massTransfer;
                    eqn.boundaries()(col) += massTransfer;
                }
            }

            for(const BoundaryLink &bd: cell.boundaries())
            {
                Polygon fluxPgn = computeFluxPolygon(bd, u.faces()[bd.face().id()], timeStep, componentNo);

                if(fluxPgn.isEmpty())
                    continue;

               Scalar massTransfer;

                switch(field.boundaryType(bd.face().id()))
                {
                case ScalarFiniteVolumeField::FIXED:
                    massTransfer = fluxPgn.area()*field.faces()[bd.face().id()];

                    break;

                case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                    if(!plicPgn.isEmpty())
                        massTransfer = intersectionPolygon(plicPgn, fluxPgn).area(); // Flux depends on the plic reconstruction in the boundary cell
                    else
                        massTransfer = fluxPgn.area()*field[cell.id()];

                    break;

                default:
                    throw Exception("plic", "div", "unrecognized or unspecified boundary type.");
                }

                eqn.boundaries()(row) -= massTransfer;
            }
        }

        for(const Cell &cell: field.grid.fluidCells())
        {
            const size_t row = cell.globalIndex();

            Scalar centralCoeff = cell.volume();
            //eqn.boundaries()(row) += std::min(std::max(field[cell.id()], 0.), 1.)*cell.volume(); // A limiting operation to ensure valid masses
            eqn.boundaries()(row) += field[cell.id()]*cell.volume();

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, centralCoeff));
        }

        eqn.matrix().assemble(entries);

        if(componentNo == 0)
            eqn.solve(); // Solve and compute the next component
        else
            return eqn;
    }
}

Polygon computeInterfacePolygon(const Cell &cell, Scalar& gamma,  const Vector2D& interfaceNormal)
{
    const int maxIters = 100;
    Polygon plicPgn;

    if(gamma >= 1)
    {
        return cell.shape();
    }
    else if (gamma <= 1e-12)
    {
        return Polygon();
    }
    else if (interfaceNormal == Vector2D(0., 0.)) // This is bad and should not happen, interface cells should have a gradient. Plic pgn will be placed in the middle4
    {
        Scalar a = 0., b = 1.;

        int iters = 0;
        while(true)
        {
            Scalar c = (a + b)/2.;

            plicPgn = cell.shape();
            plicPgn.scale(c);

            Scalar f = plicPgn.area()/cell.volume() - gamma;

            if(fabs(f) < 1e-12)
                return plicPgn;
            else if(++iters > maxIters)
                break;
            else if(f > 0)
                b = c;
            else
                a = c;
        }
    }
    else
    {
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

            cMax += 1.01*fabs(cMax);
            cMin -= 1.01*fabs(cMin);
        };

        Scalar func, c, cMin, cMax;
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

        plicPgn = clipPolygon(cell.shape(), plic.adjust(c)); // Return the clipped polygon that represents the shape of phase B
    }

    if(!plicPgn.isValid())
    {
        return Polygon();
        //throw Exception("plic", "computeInterfacePolygon", "an invalid polygon was detected, gamma = " + std::to_string(gamma));
    }

    if(fabs(plicPgn.area()/cell.volume() - gamma) > 1e-5)
        throw Exception("plic", "computeInterfacePolygon", "an invalid PLIC polygon reconstruction was detected. Reconstructed gamma = "
                        + std::to_string(plicPgn.area()/cell.volume()) + ", actual gamma = " + std::to_string(gamma) + ".");

    return plicPgn;
}

Polygon computeFluxPolygon(const BoundaryLink &link, const Vector2D& uf, Scalar timeStep, int componentNo)
{
    Vector2D off = dot(uf, link.outwardNorm())*link.outwardNorm()/link.outwardNorm().magSqr()*timeStep;
    off = componentNo == 0 ? Vector2D(off.x, 0.) : Vector2D(0., off.y);

    if(off.magSqr() < std::numeric_limits<Scalar>::epsilon())
        return Polygon();

    std::vector<Point2D> verts = {
        link.face().lNode(),
        link.face().lNode() - off,
        link.face().rNode() - off,
        link.face().rNode(),
    };

    Polygon fluxPgn = Polygon(verts);

    if(!fluxPgn.isValid())
        throw Exception("plic", "computeFluxPolygon", "an invalid polygon was detected.");

    return fluxPgn;
}
}
