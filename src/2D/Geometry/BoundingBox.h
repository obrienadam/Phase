#ifndef PHASE_BOUNDING_BOX_H
#define PHASE_BOUNDING_BOX_H

#include <vector>

#include "Point2D.h"

class BoundingBox {
public:
  enum Quadrant { I, II, III, IV };

  BoundingBox() {}

  BoundingBox(const Point2D &lBound, const Point2D &uBound);

  template <class const_iterator>
  BoundingBox(const_iterator begin, const_iterator end) {
    if (begin != end)
      lBound_ = uBound_ = *begin;

    for (auto itr = begin; itr != end; ++itr) {
      if (itr->x < lBound_.x)
        lBound_.x = itr->x;

      if (itr->y < lBound_.y)
        lBound_.y = itr->y;

      if (itr->x > uBound_.x)
        uBound_.x = itr->x;

      if (itr->y > uBound_.y)
        uBound_.y = itr->y;
    }

    center_ = (lBound_ + uBound_) / 2.;
  }

  Quadrant getQuadrant(const Point2D &point) const;

  bool isInBox(const Point2D &point) const;

  const Point2D &lBound() const { return lBound_; }

  const Point2D &uBound() const { return uBound_; }

  const Point2D &center() const { return center_; }

private:
  Point2D center_, lBound_, uBound_;
};

#endif
