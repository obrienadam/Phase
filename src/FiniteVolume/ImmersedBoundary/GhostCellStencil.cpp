#include "GhostCellStencil.h"
#include "FiniteVolumeGrid2D.h"

GhostCellStencil::GhostCellStencil(const Cell &cell, const Shape2D &shape, const FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryStencil(cell)
{
    bp_ = shape.nearestIntersect(cell_.centroid());
    ip_ = 2. * bp_ - cell.centroid();
    ipCells_ = grid.findNearestNode(ip_).cells();

    Point2D x1 = ipCells_[0].get().centroid();
    Point2D x2 = ipCells_[1].get().centroid();
    Point2D x3 = ipCells_[2].get().centroid();
    Point2D x4 = ipCells_[3].get().centroid();

    Matrix A = inverse(Matrix(4, 4, {
            x1.x*x1.y, x1.x, x1.y, 1.,
            x2.x*x2.y, x2.x, x2.y, 1.,
            x3.x*x3.y, x3.x, x3.y, 1.,
            x4.x*x4.y, x4.x, x4.y, 1.,
    }));

    Matrix m = Matrix(1, 4, {ip_.x*ip_.y, ip_.x, ip_.y, 1.}) * A;

    ipCoeffs_ = (Matrix(1, 4, {ip_.x*ip_.y, ip_.x, ip_.y, 1.}) * A).containerCopy();
}

Scalar GhostCellStencil::ipValue(const ScalarFiniteVolumeField &field) const
{
    std::vector<Scalar> vals;

    for (const Cell &cell: ipCells_)
        vals.push_back(field(cell));

    return (*interpolator_)(vals, ip_);
}

Vector2D GhostCellStencil::ipValue(const VectorFiniteVolumeField &field) const
{
    std::vector<Vector2D> vals;

    for (const Cell &cell: ipCells_)
        vals.push_back(field(cell));

    return (*interpolator_)(vals, ip_);
}

Vector2D GhostCellStencil::ipGrad(const ScalarFiniteVolumeField &field) const
{
    Matrix mat(4, 3), b(4, 1);

    int i = 0;
    Scalar ipVal = ipValue(field);
    for (const Cell &cell: ipCells_)
    {
        mat(i, 0) = cell.centroid().x - ip_.x;
        mat(i, 1) = cell.centroid().y - ip_.y;
        b(i, 0) = field(cell) - ipVal;
        ++i;
    }

    mat.solve(b);
    return Vector2D(b(0, 0), b(0, 1));
}
