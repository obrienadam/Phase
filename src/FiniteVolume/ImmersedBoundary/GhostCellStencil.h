#ifndef GHOST_CELL_STENCIL_H
#define GHOST_CELL_STENCIL_H

#include "Cell.h"
#include "Interpolation.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class GhostCellStencil
{
public:

    GhostCellStencil(const Cell &cell, const Shape2D &shape, const FiniteVolumeGrid2D &grid);

    const Cell &cell() const
    { return cell_; }

    const Point2D &imagePoint() const
    { return ip_; }

    const Point2D &boundaryPoint() const
    { return bp_; }

    Scalar length() const
    { return (ip_ - cell_.centroid()).mag(); }

    const std::vector<Ref<const Cell>> &ipCells() const
    { return ipCells_; }

    Scalar ipValue(const ScalarFiniteVolumeField &field) const;

    Vector2D ipValue(const VectorFiniteVolumeField &field) const;

    const std::vector<Scalar> &ipCoeffs() const
    { return ipCoeffs_; }

private:

    const Cell &cell_;

    Point2D ip_, bp_;
    std::vector<Ref<const Cell>> ipCells_;
    std::unique_ptr<Interpolation> interpolator_;
    std::vector<Scalar> ipCoeffs_;
};

#endif
