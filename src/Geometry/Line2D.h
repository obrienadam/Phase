#ifndef LINE_2D
#define LINE_2D

#include "Vector2D.h"
#include "Point2D.h"

class Line2D
{
public:

    static Point2D intersection(const Line2D& lineA, const Line2D& lineB);

    Line2D() {}
    Line2D(const Point2D& r0, const Point2D& d);

    Point2D operator()(Scalar t) const;

    bool isApproximatelyOnLine(const Point2D& pt) const;
    bool isAboveLine(const Point2D& pt) const;
    bool isBelowLine(const Point2D& pt) const;
    Scalar t(const Point2D& pt) const;

private:

    Point2D r0_;
    Vector2D d_, n_;

    friend std::ostream& operator<<(std::ostream& os, const Line2D& line);
};

std::ostream& operator<<(std::ostream& os, const Line2D& line);

#endif
