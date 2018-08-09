#ifndef PHASE_FACE_H
#define PHASE_FACE_H

#include "3D/Geometry/Rectangle.h"

class StructuredGrid3D;

class Face
{
public:

    enum Direction{I_POS, J_POS, K_POS, I_NEG, J_NEG, K_NEG};

    enum Index{I, J, K};

    static Index idx(Direction dir);

    static const std::array<Direction, 6> DIRECTIONS;

    Face(const StructuredGrid3D &grid, Direction f, Label i, Label j, Label k);

    const Rectangle &shape() const
    { return _shape; }

    Index idx() const
    { return _idx; }

    Label id() const
    { return _id; }

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

    bool isBoundaryFace() const
    { return _isBoundaryFace; }

protected:

    const StructuredGrid3D &_grid;

    Index _idx;

    Label _id, _i, _j, _k;

    bool _isBoundaryFace;

    Rectangle _shape;
};

#endif
