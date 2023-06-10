#ifndef PHASE_INTERSECTIONS_H
#define PHASE_INTERSECTIONS_H

#include "LineSegment3D.h"
#include "Plane.h"
#include "Ray3D.h"
#include "Sphere.h"

std::vector<Point3D> intersections(const Sphere &sphere,
                                   const LineSegment3D &ln);

std::vector<Point3D> intersections(const Sphere &sphere, const Ray3D &ray);

std::pair<Point3D, bool> intersection(const Plane &plane, const Ray3D &ray);

std::pair<Point3D, bool> intersection(const Plane &plane,
                                      const LineSegment3D &ln);

#endif
