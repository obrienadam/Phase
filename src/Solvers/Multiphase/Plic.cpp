#include "Plic.h"
#include "BisectionSearch.h"
#include "SecantSearch.h"

namespace plic
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, const VectorFiniteVolumeField& gradField, ScalarFiniteVolumeField &field, Scalar timeStep, std::vector<Polygon> &plicPolygons)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;

    for(int componentNo = 0; componentNo < 2; ++componentNo)
    {
        Equation<ScalarFiniteVolumeField> eqn(field);
        entries.clear();
        entries.reserve(4*field.grid.nActiveCells());

        for(const Cell &cell: field.grid.fluidCells())
        {
            const size_t row = cell.globalIndex();

            Polygon plicPgn = computeInterfacePolygon(cell, field[cell.id()], -gradField[cell.id()]);
            if(componentNo == 0)
                plicPolygons[cell.id()] = plicPgn;

            for(const InteriorLink &nb: cell.neighbours())
            {
                const size_t col = nb.cell().globalIndex();

                if(dot(u.faces()[nb.face().id()], nb.outwardNorm()) > 0.)
                {
                    Polygon fluxPgn = computeFluxPolygon(nb, u.faces()[nb.face().id()], timeStep, componentNo);

                    if(fluxPgn.isEmpty())
                        continue;

                    if(!plicPgn.isEmpty())
                    {
                        Scalar massTransfer = intersectionPolygon(plicPgn, fluxPgn).area();
                        eqn.boundaries()(row) -= massTransfer;
                        eqn.boundaries()(col) += massTransfer;
                    }
                    else
                    {
                        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, fluxPgn.area()));
                        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(col, row, -fluxPgn.area()));
                    }
                }
            }

            for(const BoundaryLink &bd: cell.boundaries())
            {
                Polygon fluxPgn = computeFluxPolygon(bd, u(bd.face()), timeStep, componentNo);

                if(fluxPgn.isEmpty())
                    continue;

                Scalar massTransfer = 0.;

                switch(field.boundaryType(bd.face()))
                {
                case ScalarFiniteVolumeField::FIXED:
                    massTransfer = fluxPgn.area()*field(bd.face());

                    break;

                case ScalarFiniteVolumeField::NORMAL_GRADIENT:
                    if(!plicPgn.isEmpty())
                        massTransfer = intersectionPolygon(plicPgn, fluxPgn).area(); // Flux depends on the plic reconstruction in the boundary cell
                    else
                        entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, row, fluxPgn.area()));

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

        eqn.assemble(entries);

        if(componentNo == 0)
            eqn.solve(); // Solve and compute the next component
        else
            return eqn;
    }

    return Equation<ScalarFiniteVolumeField>(field);
}

Polygon computeInterfacePolygon(const Cell &cell, Scalar& gamma,  const Vector2D& interfaceNormal)
{
    Polygon plicPgn;

    if(gamma >= 1)
    {
        return cell.shape();
    }
    else if (gamma <= 1e-12)
    {
        return Polygon();
    }
    else if (interfaceNormal == Vector2D(0., 0.)) // Invalid normal, must use a different method
    {
        std::function<Scalar(Scalar)> f = [&cell, &gamma](Scalar c){ return cell.shape().scale(c).area()/cell.volume() - gamma; };
        auto c = bisectionSearch(f, 0., 1., 1e-12, 100);

        plicPgn = cell.shape().scale(c.first); // Return a scaled concentric polygon
    }
    else
    {
        //- Create a plic line
        Line2D plic(cell.centroid(), interfaceNormal);

        //- Locate the bounds
        Scalar cMin, cMax;

        bool init = true;
        for(const Point2D& vtx: cell.shape().vertices())
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

        std::function<Scalar(Scalar)> f = [&cell, &plic, &gamma](Scalar c){ return clipPolygon(cell.shape(), plic.adjust(c)).area()/cell.volume() - gamma; };
        auto c = bisectionSearch(f, cMin, cMax, 1e-12, 100);

        plicPgn = clipPolygon(cell.shape(), plic.adjust(c.first)); // Return the clipped polygon that represents the shape of phase B
    }

    if(!plicPgn.isValid())
    {
        return Polygon();
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
