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
    //- find the bounds
    bracket(cMin, cMax);

    //- bisection
//    while(true)
//    {
//        c = (cMin + cMax)/2;
//        func = f(c);

//        if (fabs(func) < 1e-10)
//            break;
//        else if (func < 0)
//            cMin = c;
//        else
//            cMax = c;
//    }

    //- Newton
    Scalar c1, c0, func0;

    //- Initialize using a bisection iteration
    c = (cMin + cMax)/2;
    func = f(c);

    if(func < 0)
        cMin = c;
    else
        cMax = c;

    c0 = c;
    c = (cMin + cMax)/2;
    func0 = func;
    func = f(c);

    //- Start the Newton loop
    while(true)
    {
        c1 = c - func/((func - func0)/(c - c0));

        c0 = c;
        c = c1;
        func0 = func;
        func = f(c);

        if(fabs(func) < 1e-14)
            break;
    }

    return clipPolygon(cell.shape(), plic.adjust(c)); // Return the clipped polygon that represents the shape of phase B
}

}
