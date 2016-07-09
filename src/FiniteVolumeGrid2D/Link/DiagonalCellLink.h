#ifndef DIAGONAL_CELL_LINK_H
#define DIAGONAL_CELL_LINK_H

#include "Link.h"

class DiagonalCellLink : public Link
{
    DiagonalCellLink(const Cell& self, const Cell& cell);

    const Cell& cell() const { return cell_; }
    const Vector2D& rCellVec() const { return rCellVec_; }

protected:

    Ref<const Cell> cell_;
    Vector2D rCellVec_;
};

#endif
