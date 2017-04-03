#ifndef BOUNDARY_FACE_LINK_H
#define BOUNDARY_FACE_LINK_H

#include "Link.h"

class Face;

class BoundaryLink : public Link
{
public:

    BoundaryLink(const Cell &self, const Face &face);

    explicit BoundaryLink(const BoundaryLink &other);

    virtual void init();

    const Face &face() const;

    const Vector2D &rFaceVec() const
    { return rFaceVec_; }

    const Vector2D &outwardNorm() const
    { return outwardNorm_; }

protected:

    const Face &face_;
    Vector2D rFaceVec_, outwardNorm_;
};

#endif
