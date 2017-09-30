#ifndef CUT_CELL_LINK_H
#define CUT_CELL_LINK_H

#include "InteriorFaceLink.h"
#include "Shape2D.h"
#include "Cell.h"

class CutCellLink : public InteriorLink
{
public:

    CutCellLink(const Cell &self, const Face &face, const Cell &cell);

    CutCellLink(const Cell &self, const Face &face, const Cell &cell, const Shape2D &shape);

    Scalar zeta() const
    { return zeta_; }

    const Vector2D& xc() const
    { return xc_; }

    Vector2D fluidNorm() const
    { return zeta_ * outwardNorm_; }

    Vector2D solidNorm() const
    { return (1. - zeta_) * outwardNorm_; }

    const Vector2D &fluidRFaceVec() const
    { return fluidRFaceVec_; }

    const Vector2D &solidRFaceVec() const
    { return solidRFaceVec_; }

    Vector2D solidFaceCentroid() const
    { return solidFace_.center(); }

    Vector2D fluidFaceCentroid() const
    { return fluidFace_.center(); }

private:

    Scalar zeta_;
    Point2D xc_;
    LineSegment2D fluidFace_, solidFace_;
    Vector2D fluidRFaceVec_, solidRFaceVec_;
};


#endif
