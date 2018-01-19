#ifndef IMMERSED_BOUNDARY_STENCIL_H
#define IMMERSED_BOUNDARY_STENCIL_H

#include "Cell.h"
#include "ImmersedBoundaryObject.h"

class ImmersedBoundaryStencil
{
public:

    ImmersedBoundaryStencil(const Cell &cell) : cell_(cell)
    {}

    const Cell &cell() const
    { return cell_; }

    std::vector<Ref<const Cell>>::iterator begin()
    { return cells_.begin(); }

    std::vector<Ref<const Cell>>::iterator end()
    { return cells_.end(); }

    std::vector<Ref<const Cell>>::const_iterator begin() const
    { return cells_.begin(); }

    std::vector<Ref<const Cell>>::const_iterator end() const
    { return cells_.end(); }

    const std::vector<Scalar> &dirichletCoeffs() const
    { return dirichletCoeffs_; }

    const std::vector<Scalar> &neumannCoeffs() const
    { return neumannCoeffs_; }

    const std::vector<Point2D> &boundaryPoints() const
    { return boundaryPoints_; }

    const std::vector<Scalar> &dirichletBoundaryCoeffs() const
    { return dirichletBoundaryCoeffs_; }

    const std::vector<Scalar> &neumannBoundaryCoeffs() const
    { return neumannBoundaryCoeffs_; }

protected:

    Ref<const Cell> cell_;
    std::vector<Ref<const Cell>> cells_;
    std::vector<Scalar> dirichletCoeffs_;
    std::vector<Scalar> neumannCoeffs_;

    std::vector<Point2D> boundaryPoints_;
    std::vector<Scalar> dirichletBoundaryCoeffs_;
    std::vector<Scalar> neumannBoundaryCoeffs_;
};

#endif
