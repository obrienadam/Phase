#ifndef PHASE_INTERSECTION_H
#define PHASE_INTERSECTION_H

#include "System/StaticVector.h"

#include "Circle.h"
#include "Line2D.h"
#include "LineSegment2D.h"
#include "Ray2D.h"

std::pair<Point2D, bool> intersection(const Line2D &lnA, const Line2D &lnB);

std::pair<Point2D, bool> intersection(const Ray2D &ray,
                                      const LineSegment2D &ln);

std::pair<LineSegment2D, bool> intersection(const Line2D &lnA, const Circle &c);

StaticVector<Point2D, 2> intersection(const Circle &c, const Ray2D &r);

#endif
