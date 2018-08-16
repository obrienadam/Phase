#ifndef PHASE_CELL_LINK_H
#define PHASE_CELL_LINK_H

#include "Geometry/Point2D.h"

#include "Link.h"

template<class T>
class FiniteVolumeField;

class CellLink : public Link
{
public:

    CellLink(const Cell &self, const Cell &other);

    const Cell &cell() const
    { return cell_; }

    const Vector2D &rCellVec() const
    { return rCellVec_; }

    Scalar alpha(const Point2D &pt) const;

    Scalar linearInterpolate(const FiniteVolumeField<Scalar> &phi, const Point2D &pt) const;

    Vector2D linearInterpolate(const FiniteVolumeField<Vector2D> &u, const Point2D &pt) const;

protected:

    const Cell &cell_;
    Vector2D rCellVec_;
};


#endif
