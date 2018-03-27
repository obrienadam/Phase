#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <vector>

#include "Point2D.h"

class BoundingBox
{
public:
    enum Quadrant{I, II, III, IV};

    BoundingBox() {}
    BoundingBox(const Point2D& lBound, const Point2D& uBound);
    BoundingBox(const Point2D* points, size_t nPoints);

    Quadrant getQuadrant(const Point2D &point) const;
    bool isInBox(const Point2D& point) const;

    const Point2D& lBound() const { return lBound_; }
    const Point2D& uBound() const { return uBound_; }
    const Point2D& center() const { return center_; }

    std::string toString() const;

private:

    Point2D center_, lBound_, uBound_;
};

#endif
