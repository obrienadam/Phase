#ifndef CELL_LINK_H
#define CELL_LINK_H

#include "Link.h"

class CellLink : public Link
{
public:

    CellLink(const Cell &self, const Cell &other);

    const Cell &cell() const
    { return cell_; }

    const Vector2D &rCellVec() const
    { return rCellVec_; }

protected:

    const Cell &cell_;
    Vector2D rCellVec_;
};


#endif