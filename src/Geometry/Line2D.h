#ifndef LINE_2D
#define LINE_2D

#include "Vector2D.h"
#include "Point2D.h"

class Line2D
{
public:

    Line2D(const Vector2D& normal, Scalar c) : n_(normal.unitVec()), c_(c) {}

private:

    Vector2D n_;
    Scalar c_;
};

#endif
