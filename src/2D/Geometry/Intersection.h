#ifndef PHASE_INTERSECTION_H
#define PHASE_INTERSECTION_H

#include <array>

#include "Point2D.h"
#include "Line2D.h"
#include "Ray2D.h"
#include "LineSegment2D.h"
#include "Circle.h"

template<Size Nm>
class Intersection
{
public:

    Intersection() {}

    Intersection(const std::initializer_list<Point2D> &intersections)
    {
        std::copy(intersections.end(), intersections.end(), intersections_.begin());
        nIntersections_ = intersections.size();
    }

    void add(const Point2D &intersection)
    { intersections_[nIntersections_++] = intersection; }

    bool empty() const
    { return nIntersections_ == 0; }

    typename std::array<Point2D, Nm>::const_iterator begin() const
    { return intersections_.begin(); }

    typename std::array<Point2D, Nm>::const_iterator end() const
    { return begin() + nIntersections_; }

protected:

    std::array<Point2D, Nm> intersections_;

    Size nIntersections_ = 0;
};

std::pair<Point2D, bool> intersection(const Line2D &lnA, const Line2D &lnB);

std::pair<Point2D, bool> intersection(const Ray2D &ray, const LineSegment2D &ln);

std::pair<LineSegment2D, bool> intersection(const Line2D &lnA, const Circle &c);

Intersection<2> intersection(const Ray2D &ray, const Circle &circle);

#endif
