#ifndef INTERIOR_FACE_LINK_H
#define INTERIOR_FACE_LINK_H

#include "BoundaryFaceLink.h"

class Face;

class InteriorLink : public BoundaryLink
{
public:

    InteriorLink(const Cell& self, const Face& face, const Cell& cell);
    explicit InteriorLink(const InteriorLink& other);

    InteriorLink& operator=(const InteriorLink& rhs);

    const Cell& cell() const;

    const Vector2D& rCellVec() const { return rCellVec_; }

protected:

    const Cell& cell_;
    Vector2D rCellVec_;
};

#endif
