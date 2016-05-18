#include "Plic.h"

namespace plic
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field)
{
    std::vector<Equation<ScalarFiniteVolumeField>::Triplet> entries;
    Equation<ScalarFiniteVolumeField> eqn(field);
    VectorFiniteVolumeField gradField = grad(field); //- should this be smoothed?

    entries.reserve(5*field.grid.nActiveCells());

    const Scalar toler = 1e-1;
    // Figure out the plic reconstruction
    for(const Cell &cell: field.grid.fluidCells())
    {
        bool isInterfaceCell = field[cell.id()] > toler && field[cell.id()] < 1. - 1e-10;

        Scalar centralCoeff = 0.;
        const size_t row = cell.globalIndex();

        for(const InteriorLink &nb: cell.neighbours())
        {
            size_t col = nb.cell().globalIndex();
            Scalar coeff;

            if(isInterfaceCell)
            {
                Vector2D norm = -gradField[cell.id()];
                Line2D line(cell.centroid(), norm.normalVec());
            }
            else
            {
                Scalar faceFlux = dot(u.faces()[nb.face().id()], nb.outwardNorm());

                coeff = std::min(faceFlux, 0.);
                centralCoeff += std::max(faceFlux, 0.);
            }

            entries.push_back(Equation<ScalarFiniteVolumeField>::Triplet(row, col, coeff));
        }

        for(const BoundaryLink &bd: cell.boundaries())
        {

        }
    }

    eqn.matrix().assemble(entries);
    return eqn;
}

Polygon computeInterfacePolygon(const Cell &cell, Scalar gamma,  const Vector2D& interfaceNormal)
{
    Line2D plic(cell.centroid(), interfaceNormal);

    auto f = [&cell, &plic, &gamma](Scalar c)
    {
        return clipPolygon(cell.shape(), Line2D(plic.r0() + c*plic.n(), plic.n())).area()/cell.volume() - gamma;
    };

    auto bracket = [&cell, &plic](Scalar &minC, Scalar &maxC)
    {
        minC = maxC = 0;

        for(const Point2D &vtx: cell.shape())
        {
            Scalar c = dot(vtx - plic.r0(), plic.n());

            maxC = c > maxC ? c : maxC;
            minC = c < minC ? c : minC;
        }
    };
    Scalar minC, maxC, c0, f0;
    bool converged;

    //- Start newton
    bracket(minC, maxC);
    converged = true;

    //- Do one iteration of bisection to get an initial estimate
    c0 = (minC + maxC)/2;
    f0 = f(c0);

    if(f0 < 0)
        minC = c0;
    else
        maxC = c0;

    //- Init the vals for the Newton method
    Scalar c00 = c0;
    c0 = (minC + maxC)/2;
    Scalar f00 = f0;
    f0 = f(c0);
    //- With initial guesses, begin a Newton root finding tecnique with a numerical approximation to the derivative
    const int maxIters = 100;
    int i = 0;
    do
    {
        Scalar c1 = c0 - f0/((f0 - f00)/(c0 - c00));

        c00 = c0;
        c0 = c1;
        f00 = f0;
        f0 = f(c0);

        if(++i == maxIters)
        {
            converged = false;
            break;
        }

    } while(fabs(f0 - f00) > 1e-14 && fabs(c0 - c00) > 1e-14 );


    //- Try golden section search
    if(!converged)
    {
        converged = true;
        bracket(minC, maxC);

        const Scalar gr = (sqrt(5.) - 1.)/2.;
        Scalar fc, fd, c, d;
        i = 0;
        do
        {
            c = maxC - gr*(maxC - minC);
            d = minC + gr*(maxC - minC);

            //- value of c0 that minimizes this absolute error function is the solution!
            fc = fabs(f(c));
            fd = fabs(f(d));

            if (fc < fd)
                maxC = d;
            else
                minC = c;
                c0 = d;

            if(++i == maxIters)
            {
                converged = false;
                break;
            }

        } while(std::min(fc, fd) > 1e-10);

        c0 = fc < fd ? c : d;
        f0 = f(c0);

    } // End golden section search

    //- Try bisection search

    if(!converged)
    {
        converged = true;
        bracket(minC, maxC);
        i = 0;
        do
        {
            c0 = (minC + maxC)/2;

            f0 = f(c0);

            if(f0 > 0)
                maxC = c0;
            else
                minC = c0;

            if(++i == maxIters)
            {
                converged = false;
                break;
            }

        } while (fabs(f0) > 1e-10);
    } // end bisection search

    return clipPolygon(cell.shape(), Line2D(plic.r0() + c0*plic.n(), plic.n())); // Return the clipped polygon that represents the shape of phase B
}

}
