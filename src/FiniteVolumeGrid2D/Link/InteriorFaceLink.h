#ifndef INTERIOR_FACE_LINK_H
#define INTERIOR_FACE_LINK_H

#include "BoundaryFaceLink.h"

class InteriorLink : public BoundaryLink
{
public:

    enum Direction{INGOING, OUTGOING};

    InteriorLink(const Cell& self, const Face& face, const Cell& cell, Direction direction);

    const Cell& cell() const { return cell_; }
    const Vector2D& rCellVec() const { return rCellVec_; }

    Direction direction() const { return direction_; }

protected:

    Ref<const Cell> cell_;
    Vector2D rCellVec_;
    Direction direction_;
};

#endif
