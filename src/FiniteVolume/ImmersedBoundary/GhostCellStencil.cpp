#include "GhostCellStencil.h"
#include "GhostCellImmersedBoundaryObject.h"

GhostCellStencil::GhostCellStencil(const Cell &cell,
                                   const GhostCellImmersedBoundaryObject &ibObj,
                                   const FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryStencil(cell)
{
    bp_ = ibObj.shape().nearestIntersect(cell_.centroid());
    ip_ = 2. * bp_ - cell.centroid();
    cells_ = grid.findNearestNode(ip_).cells();

    if (cells_.size() != 4)
        throw Exception("GhostCellStencil", "GhostCellStencil", "number of image point cells must be 4.");

    Point2D x1 = cells_[0].get().centroid();
    Point2D x2 = cells_[1].get().centroid();
    Point2D x3 = cells_[2].get().centroid();
    Point2D x4 = cells_[3].get().centroid();

    bool ghostCellInStencil = false;
    for (const Cell &cell: cells_)
        if (cell.id() == cell_.id())
        {
            ghostCellInStencil = true;
            break;
        }

    A_ = Matrix(4, 4, {
            x1.x * x1.y, x1.x, x1.y, 1.,
            x2.x * x2.y, x2.x, x2.y, 1.,
            x3.x * x3.y, x3.x, x3.y, 1.,
            x4.x * x4.y, x4.x, x4.y, 1.,
    }).invert();

    if (ghostCellInStencil)
    {
        Vector2D n = ibObj.nearestEdgeNormal(bp_).unitVec();
        dirichletCoeffs_ = Matrix(1, 4, {bp_.x * bp_.y, bp_.x, bp_.y, 1.}) * A_;
        neumannCoeffs_ = Matrix(1, 4, {bp_.y * n.x + bp_.x * n.y, n.x, n.y, 0.}) * A_;
    }
    else
    {
        cells_.push_back(std::cref(cell_));

        dirichletCoeffs_ = 0.5 * Matrix(1, 4, {ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * A_;
        dirichletCoeffs_.push_back(0.5);

        neumannCoeffs_ = Matrix(1, 4, {ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * A_ / length();
        neumannCoeffs_.push_back(-1. / length());
    }
}

GhostCellStencil::GhostCellStencil(const Cell &cell,
                                   const GhostCellImmersedBoundaryObject &ibObj,
                                   const FiniteVolumeGrid2D &grid,
                                   Scalar theta,
                                   const Vector2D &m)
        :
        ImmersedBoundaryStencil(cell)
{
    bp_ = ibObj.shape().nearestIntersect(cell_.centroid());
    ip_ = 2. * bp_ - cell.centroid();
    cells_ = grid.findNearestNode(ip_).cells();

    if (cells_.size() != 4)
        throw Exception("GhostCellStencil", "GhostCellStencil", "number of image point cells must be 4.");

    Point2D x1 = cells_[0].get().centroid();
    Point2D x2 = cells_[1].get().centroid();
    Point2D x3 = cells_[2].get().centroid();
    Point2D x4 = cells_[3].get().centroid();

    bool ghostCellInStencil = false;
    for (const Cell &cell: cells_)
        if (cell.id() == cell_.id())
        {
            ghostCellInStencil = true;
            break;
        }

    A_ = Matrix(4, 4, {
            x1.x * x1.y, x1.x, x1.y, 1.,
            x2.x * x2.y, x2.x, x2.y, 1.,
            x3.x * x3.y, x3.x, x3.y, 1.,
            x4.x * x4.y, x4.x, x4.y, 1.,
    }).invert();

    Vector2D n = ibObj.nearestEdgeNormal(bp_).unitVec();
    Vector2D t = n.tangentVec().unitVec();
    Vector2D cl = t.rotate(dot(m, t) > 0. ? -theta: theta);

    if (ghostCellInStencil)
    {
        neumannCoeffs_ = Matrix(1, 4, {bp_.y * cl.x + bp_.x * cl.y, cl.x, cl.y, 0.}) * A_;
    }
    else
    {
        cells_.push_back(std::cref(cell_));

        Vector2D ip2 = bp_ + (ip_ - bp_).rotate(dot(t, m) > 0 ? M_PI_2 - theta: theta - M_PI_2);
        Scalar l = (ip_ - bp_).mag() / 2.;

        auto extraCells = grid.findNearestNode(ip2).cells();
        cells_.insert(cells_.end(), extraCells.begin(), extraCells.end());

        Point2D x21 = extraCells[0].get().centroid();
        Point2D x22 = extraCells[1].get().centroid();
        Point2D x23 = extraCells[2].get().centroid();
        Point2D x24 = extraCells[3].get().centroid();

        Matrix A2 = Matrix(4, 4, {
                x21.x*x21.y, x21.x, x21.y, 1.,
                x22.x*x22.y, x22.x, x22.y, 1.,
                x23.x*x23.y, x23.x, x23.y, 1.,
                x24.x*x24.y, x24.x, x24.y, 1.,
        }).invert();

        auto coeffs1 = -0.5 * Matrix(1, 4, {ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * A_;
        auto coeffs2 = Matrix(1, 4, {ip2.x * ip2.y, ip2.x, ip2.y, 1.}) * A2;

        //- Boundary point equation
        neumannCoeffs_.insert(neumannCoeffs_.begin(), coeffs1.begin(), coeffs1.end());
        neumannCoeffs_.push_back(-0.5);

        //- Extrapolation point formula
        neumannCoeffs_.insert(neumannCoeffs_.end(), coeffs2.begin(), coeffs2.end());
    }
}

Scalar GhostCellStencil::ipValue(const ScalarFiniteVolumeField &field) const
{
    Matrix c = A_ * Matrix(4, 1, {
            field(cells_[0]),
            field(cells_[1]),
            field(cells_[2]),
            field(cells_[3])
    });

    return c(0, 0) * ip_.x * ip_.y + c(1, 0) * ip_.x + c(2, 0) * ip_.y + c(3, 0);
}

Vector2D GhostCellStencil::ipValue(const VectorFiniteVolumeField &field) const
{
    Matrix c = A_ * Matrix(4, 2, {
            field(cells_[0]).x, field(cells_[0]).y,
            field(cells_[1]).x, field(cells_[1]).y,
            field(cells_[2]).x, field(cells_[2]).y,
            field(cells_[3]).x, field(cells_[3]).y
    });

    return Vector2D(
            c(0, 0) * ip_.x * ip_.y + c(1, 0) * ip_.x + c(2, 0) * ip_.y + c(3, 0),
            c(0, 1) * ip_.x * ip_.y + c(1, 1) * ip_.x + c(2, 1) * ip_.y + c(3, 1)
    );
}

Vector2D GhostCellStencil::ipGrad(const ScalarFiniteVolumeField &field) const
{
    Matrix x = Matrix(2, 4, {ip_.y, 1., 0., 0., ip_.x, 0., 1., 0.}) * A_ * Matrix(4, 1, {
            field(cells_[0]),
            field(cells_[1]),
            field(cells_[2]),
            field(cells_[3])
    });

    return Vector2D(x(0, 0), x(1, 0));
}
