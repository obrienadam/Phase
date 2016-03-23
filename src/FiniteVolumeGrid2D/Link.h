#ifndef LINK_H
#define LINK_H

#include "Vector2D.h"
#include "Types.h"

class Cell;
class Face;

class Link
{
public:

    Link(const Cell& self) : self_(self) {}

    const Cell& self() const { return self_; }

protected:

    Ref<const Cell> self_;
};

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

class InteriorLink : public BoundaryLink
{
public:

    InteriorLink(const Cell& self, const Face& face, const Cell& cell);

    const Cell& cell() const { return cell_; }
    const Vector2D& rCellVec() const { return rCellVec_; }

protected:

    Ref<const Cell> cell_;
    Vector2D rCellVec_;
};

#endif
