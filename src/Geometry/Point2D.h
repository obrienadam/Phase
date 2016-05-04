#ifndef POINT_2D_H
#define POINT_2D_H

#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry.hpp>

#include "Vector2D.h"

typedef Vector2D Point2D;

BOOST_GEOMETRY_REGISTER_POINT_2D(Point2D, Scalar, cs::cartesian, x, y)


#include <boost/geometry/geometries/polygon.hpp>

#endif
