#ifndef PHASE_POLY_LINE_2D_H
#define PHASE_POLY_LINE_2D_H

#include "Point2D.h"

class PolyLine2D {
public:
  PolyLine2D() {}

  PolyLine2D(const std::initializer_list<Point2D> &pts) : pts_(pts) {}

  template <class const_iterator>
  PolyLine2D(const_iterator begin, const_iterator end) : pts_(begin, end) {}

  //- Operators
  PolyLine2D &operator=(const std::initializer_list<Point2D> &pts);

  const Point2D &operator[](std::size_t i) const { return pts_[i]; }

  Size nVerts() const { return pts_.size(); }

  Size nSegments() const { return pts_.size() - 1; }

  Scalar length() const;

  //- Iterators
  typename std::vector<Point2D>::const_iterator begin() const {
    return pts_.begin();
  }

  typename std::vector<Point2D>::const_iterator end() const {
    return pts_.end();
  }

private:
  std::vector<Point2D> pts_;
};

#endif
