#ifndef PHASE_STATIC_POLY_LINE_2D_H
#define PHASE_STATIC_POLY_LINE_2D_H

#include "Point2D.h"

template <std::size_t N> class StaticPolyLine2D {
public:
  StaticPolyLine2D() {}

  StaticPolyLine2D(const std::initializer_list<Point2D> &pts) : pts_(pts) {}

  template <class const_iterator>
  StaticPolyLine2D(const_iterator begin, const_iterator end)
      : pts_(begin, end) {}

  //- Operators
  StaticPolyLine2D &operator=(const std::initializer_list<Point2D> &pts) {
    std::copy(pts.begin(), pts.end(), pts_.begin());
    return *this;
  }

  Point2D &operator[](std::size_t i) { return pts_[i]; }

  const Point2D &operator[](std::size_t i) const { return pts_[i]; }

  Size nVerts() const { return pts_.size(); }

  Size nSegments() const { return pts_.size() - 1; }

  Scalar length() const {
    Scalar l = 0.;

    for (std::size_t i = 1; i < N; ++i)
      l += (pts_[i] - pts_[i - 1]).mag();

    return l;
  }

  //- Iterators
  typename std::vector<Point2D>::const_iterator begin() const {
    return pts_.begin();
  }

  typename std::vector<Point2D>::const_iterator end() const {
    return pts_.end();
  }

private:
  std::array<Point2D, N> pts_;
};

#endif
