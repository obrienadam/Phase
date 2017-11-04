#ifndef INTERIOR_LINK_H
#define INTERIOR_LINK_H

#include "CellLink.h"

class Face;

class InteriorLink : public CellLink
{
public:

    InteriorLink(const Cell &self, const Face &face, const Cell &cell);

    explicit InteriorLink(const InteriorLink &other);

    InteriorLink &operator=(const InteriorLink &rhs);

    const Face &face() const
    { return face_; }

    Scalar volumeWeight() const;

    Scalar distanceWeight() const;

    Scalar distanceSqrWeight() const;

    const Vector2D &outwardNorm() const
    { return outwardNorm_; }

    const Vector2D &rFaceVec() const
    { return rFaceVec_; }

protected:

    const Face &face_;
    Vector2D outwardNorm_, rFaceVec_;
};

#endif
