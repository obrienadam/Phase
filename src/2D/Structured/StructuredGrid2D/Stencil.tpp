#include <algorithm>

#include "Stencil.h"
#include "Cell.h"

template<Size N>
Stencil<N>::Stencil(const Cell &cell, Orientation dir, const std::initializer_list<int> &offsets)
{
    std::transform(offsets.begin(), offsets.end(), _stPts.begin(), [&cell, dir](int offset)
    {
        return Point(cell.cell(dir, offset), cell.dh(dir, offset));
    });
}
