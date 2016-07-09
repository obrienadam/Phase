#ifndef BOUNDARY_FACE_LINK_H
#define BOUNDARY_FACE_LINK_H

#include "Link.h"

class BoundaryLink : public Link
{
public:

    BoundaryLink(const Cell& self, const Face& face);

    const Face& face() const { return face_; }
    const Vector2D& rFaceVec() const { return rFaceVec_; }
    const Vector2D& outwardNorm() const { return outwardNorm_; }

protected:

    Ref<const Face> face_;
    Vector2D rFaceVec_, outwardNorm_;
};

#endif
