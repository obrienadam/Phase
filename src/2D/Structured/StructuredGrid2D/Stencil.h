#ifndef PHASE_STENCIL_H
#define PHASE_STENCIL_H

#include "Orientation.h"
#include "StructuredGrid2D.h"

template<Size N>
class Stencil
{
public:

    class Point
    {
    public:

        Point() : _cell(nullptr), _dh(0.) {}

        Point(const Cell& cell, Scalar dh) : _cell(&cell), _dh(dh) {}

        const Cell &cell() const
        { return *_cell; }

        Scalar dh() const
        { return _dh; }

    private:

        const Cell *_cell;

        Scalar _dh;
    };

    typedef typename std::array<Point, N>::iterator iterator;

    typedef typename std::array<Point, N>::const_iterator const_iterator;

    Stencil(const Cell &cell, Orientation dir, const std::initializer_list<int> &offsets);

    //- iterators

    iterator begin()
    { return _stPts.begin(); }

    iterator end()
    { return _stPts.end(); }

    const_iterator cbegin() const
    { return _stPts.cbegin(); }

    const_iterator cend() const
    { return _stPts.cend(); }

protected:

    std::array<Point, N> _stPts;

};

#include "Stencil.tpp"

#endif
