#ifndef PHASE_FACE_H
#define PHASE_FACE_H

#include "3D/Geometry/Rectangle.h"

class StructuredGrid3D;

class Face
{
public:

    enum Direction{I_POS, J_POS, K_POS, I_NEG, J_NEG, K_NEG};

    Face(const StructuredGrid3D &grid, Direction f, Label i, Label j, Label k);

    const Rectangle &shape() const
    { return _shape; }

    Label i() const
    { return _i; }

    Label j() const
    { return _j; }

    Label k() const
    { return _k; }

    const Point3D &centroid() const
    { return _shape.centroid(); }

    const Vector3D &norm() const
    { return _shape.norm(); }

protected:

    const StructuredGrid3D &_grid;

    Label _id, _i, _j, _k;

    Rectangle _shape;
};

#endif
