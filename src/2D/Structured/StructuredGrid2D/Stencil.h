#ifndef PHASE_STENCIL_H
#define PHASE_STENCIL_H

#include "StructuredGrid2D.h"

class Stencil
{
public:

    class Point
    {
    public:

        Point(const Cell& cell, Scalar dh) : _cell(cell), _dh(dh) {}

    private:

        const Cell &_cell;

        Scalar _dh;
    };

    Stencil(const Cell &cell, StructuredGrid2D::CoordinateDirection dir, int fbias, int bbias);

    //- iterators

    std::vector<Point>::iterator begin()
    { return _stPts.begin(); }

    std::vector<Point>::iterator end()
    { return _stPts.end(); }

protected:

    std::vector<Point> _stPts;

};

#endif
