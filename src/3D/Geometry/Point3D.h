#ifndef PHASE_POINT_3D_H
#define PHASE_POINT_3D_H

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>

#include "Vector3D.h"

typedef Vector3D Point3D;

BOOST_GEOMETRY_REGISTER_POINT_3D(Point3D, Scalar, cs::cartesian, x, y, z)

#endif
