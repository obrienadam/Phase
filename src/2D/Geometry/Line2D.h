#ifndef PHASE_LINE_2D_H
#define PHASE_LINE_2D_H

#include "Point2D.h"
#include "Vector2D.h"

class Line2D {
public:
  Line2D();

  Line2D(const Point2D &r0, const Vector2D &n);

  Point2D operator()(Scalar t) const;

  bool isApproximatelyOnLine(const Point2D &pt) const;

  bool isAboveLine(const Point2D &pt) const;

  bool isBelowLine(const Point2D &pt) const;

  Point2D intersection(const Line2D &other) const;

  Point2D nearestPoint(const Point2D &pt) const;

  const Point2D &r0() const { return r0_; }

  const Vector2D &n() const { return n_; }

  const Vector2D &d() const { return d_; }

  const Line2D adjust(Scalar c) const { return Line2D(r0_ + c * n_, n_); }

private:
  Point2D r0_;
  Vector2D n_, d_;

  friend std::ostream &operator<<(std::ostream &os, const Line2D &line);
};

std::ostream &operator<<(std::ostream &os, const Line2D &line);

#endif
