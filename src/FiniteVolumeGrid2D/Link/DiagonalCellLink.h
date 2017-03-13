#ifndef DIAGONAL_CELL_LINK_H
#define DIAGONAL_CELL_LINK_H

#include "Link.h"

class DiagonalCellLink : public Link
{
public:

    DiagonalCellLink(const Cell &self, const Cell &cell);

    explicit DiagonalCellLink(const DiagonalCellLink &other);

    DiagonalCellLink &operator=(const DiagonalCellLink &rhs);

    const Cell &cell() const
    { return cell_; }

    const Vector2D &rCellVec() const
    { return rCellVec_; }

protected:

    const Cell &cell_;
    Vector2D rCellVec_;
};

#endif
