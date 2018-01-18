#include "GhostCellStencil.h"
#include "GhostCellImmersedBoundaryObject.h"

GhostCellStencil::GhostCellStencil(const Cell &cell,
                                   const GhostCellImmersedBoundaryObject &ibObj,
                                   const FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryStencil(cell)
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ip_ = 2. * bp_ - cell.centroid();
    nw_ = ibObj.nearestEdgeNormal(bp_).unitVec();
    cells_ = grid.findNearestNode(ip_).cells();

    if (cells_.size() != 4)
    {
        std::ostringstream sout;
        sout << "number of image point cells must be 4. Boundary point = " << bp_ << ", image point = " << ip_ << ".";
        throw Exception("GhostCellStencil", "GhostCellStencil", sout.str());
    }

    Point2D x1 = cells_[0].get().centroid();
    Point2D x2 = cells_[1].get().centroid();
    Point2D x3 = cells_[2].get().centroid();
    Point2D x4 = cells_[3].get().centroid();

    bool ghostCellInStencil = false;
    for (const Cell &cell: cells_)
        if (cell_.get().id() == cell.id())
        {
            ghostCellInStencil = true;
            break;
        }

    A_ = inverse<4, 4>(
            {
                    x1.x * x1.y, x1.x, x1.y, 1.,
                    x2.x * x2.y, x2.x, x2.y, 1.,
                    x3.x * x3.y, x3.x, x3.y, 1.,
                    x4.x * x4.y, x4.x, x4.y, 1.,
            });

    if (ghostCellInStencil)
    {
        Vector2D n = ibObj.nearestEdgeNormal(bp_).unitVec();
        auto xd = StaticMatrix<1, 4>({bp_.x * bp_.y, bp_.x, bp_.y, 1.}) * A_;
        auto xn = StaticMatrix<1, 4>({{bp_.y * n.x + bp_.x * n.y, n.x, n.y, 0.}}) * A_;

        dirichletCoeffs_.insert(dirichletCoeffs_.end(), xd.data(), xd.data() + 4);
        neumannCoeffs_.insert(neumannCoeffs_.end(), xn.data(), xn.data() + 4);
    }
    else
    {
        cells_.push_back(cell_);

        auto xd = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * A_ / 2.;
        auto xn = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * A_ / -length();

        dirichletCoeffs_.insert(dirichletCoeffs_.end(), xd.data(), xd.data() + 4);
        dirichletCoeffs_.push_back(1. / 2.);

        neumannCoeffs_.insert(neumannCoeffs_.end(), xn.data(), xn.data() + 4);
        neumannCoeffs_.push_back(1. / length());
    }
}

GhostCellStencil::GhostCellStencil(const Cell &cell,
                                   const Point2D &bp,
                                   const Vector2D &cl,
                                   const FiniteVolumeGrid2D &grid)
        :
        ImmersedBoundaryStencil(cell)
{
    bp_ = bp;
    ip_ = 2 * bp - cell.centroid();
    nw_ = (bp_ - ip_).unitVec();
    cells_ = grid.findNearestNode(ip_).cells();

    if (cells_.size() != 4)
    {
        std::ostringstream sout;
        sout << "number of image point cells must be 4. Boundary point = " << bp_ << ", image point = " << ip_ << ".";
        throw Exception("GhostCellStencil", "GhostCellStencil", sout.str());
    }

    Point2D x1 = cells_[0].get().centroid();
    Point2D x2 = cells_[1].get().centroid();
    Point2D x3 = cells_[2].get().centroid();
    Point2D x4 = cells_[3].get().centroid();

    bool ghostCellInStencil = false;
    for (const Cell &cell: cells_)
        if (cell_.get().id() == cell.id())
        {
            ghostCellInStencil = true;
            break;
        }

    A_ = inverse<4, 4>(
            {
                    x1.x * x1.y, x1.x, x1.y, 1.,
                    x2.x * x2.y, x2.x, x2.y, 1.,
                    x3.x * x3.y, x3.x, x3.y, 1.,
                    x4.x * x4.y, x4.x, x4.y, 1.,
            });

    if (ghostCellInStencil)
    {
        Vector2D n = cl.unitVec();
        auto xd = StaticMatrix<1, 4>({bp_.x * bp_.y, bp_.x, bp_.y, 1.}) * A_;
        auto xn = StaticMatrix<1, 4>({{bp_.y * n.x + bp_.x * n.y, n.x, n.y, 0.}}) * A_;

        dirichletCoeffs_.insert(dirichletCoeffs_.end(), xd.data(), xd.data() + 4);
        neumannCoeffs_.insert(neumannCoeffs_.end(), xn.data(), xn.data() + 4);
    }
    else
    {
        cells_.push_back(cell_);

        auto xd = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * A_ / 2.;
        auto xn = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * A_ / -length();

        dirichletCoeffs_.insert(dirichletCoeffs_.end(), xd.data(), xd.data() + 4);
        dirichletCoeffs_.push_back(1. / 2.);

        neumannCoeffs_.insert(neumannCoeffs_.end(), xn.data(), xn.data() + 4);
        neumannCoeffs_.push_back(1. / length());
    }
}

Scalar GhostCellStencil::ipValue(const ScalarFiniteVolumeField &field) const
{
    auto c = A_ * StaticMatrix<4, 1>(
            {
                    field(cells_[0]),
                    field(cells_[1]),
                    field(cells_[2]),
                    field(cells_[3])
            });

    return c(0, 0) * ip_.x * ip_.y + c(1, 0) * ip_.x + c(2, 0) * ip_.y + c(3, 0);
}

Vector2D GhostCellStencil::ipValue(const VectorFiniteVolumeField &field) const
{
    auto c = A_ * StaticMatrix<4, 2>(
            {
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

Scalar GhostCellStencil::bpValue(const ScalarFiniteVolumeField &field) const
{
    auto cells = field.grid().findNearestNode(bp_).cells();

    Point2D x1 = cells[0].get().centroid();
    Point2D x2 = cells[1].get().centroid();
    Point2D x3 = cells[2].get().centroid();
    Point2D x4 = cells[3].get().centroid();

    auto A = inverse<4, 4>({
                                   x1.x * x1.y, x1.x, x1.y, 1.,
                                   x2.x * x2.y, x2.x, x2.y, 1.,
                                   x3.x * x3.y, x3.x, x3.y, 1.,
                                   x4.x * x4.y, x4.x, x4.y, 1.,
                           });

    auto x = StaticMatrix<1, 4>({bp_.x * bp_.y, bp_.x, bp_.y, 1.});
    auto b = StaticMatrix<4, 1>({
                                        field(cells[0]),
                                        field(cells[1]),
                                        field(cells[2]),
                                        field(cells[3])
                                });

    return (x * A * b)(0, 0);
}

Vector2D GhostCellStencil::bpValue(const VectorFiniteVolumeField &field) const
{
    auto cells = field.grid().findNearestNode(bp_).cells();

    Point2D x1 = cells[0].get().centroid();
    Point2D x2 = cells[1].get().centroid();
    Point2D x3 = cells[2].get().centroid();
    Point2D x4 = cells[3].get().centroid();

    auto A = inverse<4, 4>({
                                   x1.x * x1.y, x1.x, x1.y, 1.,
                                   x2.x * x2.y, x2.x, x2.y, 1.,
                                   x3.x * x3.y, x3.x, x3.y, 1.,
                                   x4.x * x4.y, x4.x, x4.y, 1.,
                           });

    auto x = StaticMatrix<1, 4>({bp_.x * bp_.y, bp_.x, bp_.y, 1.});

    auto b = StaticMatrix<4, 2>({
                                        field(cells[0]).x, field(cells[0]).y,
                                        field(cells[1]).x, field(cells[1]).y,
                                        field(cells[2]).x, field(cells[2]).y,
                                        field(cells[3]).x, field(cells[3]).y
                                });

    StaticMatrix<1, 2> u = x * A * b;

    return Vector2D(u(0, 0), u(0, 1));
}

Vector2D GhostCellStencil::ipGrad(const ScalarFiniteVolumeField &field) const
{
    auto x = StaticMatrix<2, 4>({ip_.y, 1., 0., 0., ip_.x, 0., 1., 0.}) * A_ * StaticMatrix<4, 1>(
            {
                    field(cells_[0]),
                    field(cells_[1]),
                    field(cells_[2]),
                    field(cells_[3])
            });

    return Vector2D(x(0, 0), x(1, 0));
}

Vector2D GhostCellStencil::bpGrad(const ScalarFiniteVolumeField &field) const
{
    auto cells = field.grid().findNearestNode(bp_).cells();

    Point2D x1 = cells[0].get().centroid();
    Point2D x2 = cells[1].get().centroid();
    Point2D x3 = cells[2].get().centroid();
    Point2D x4 = cells[3].get().centroid();

    auto A = inverse<4, 4>(
            {
                    x1.x * x1.y, x1.x, x1.y, 1.,
                    x2.x * x2.y, x2.x, x2.y, 1.,
                    x3.x * x3.y, x3.x, x3.y, 1.,
                    x4.x * x4.y, x4.x, x4.y, 1.,
            });

    auto b = StaticMatrix<4, 1>(
            {
                    field(cells[0]),
                    field(cells[1]),
                    field(cells[2]),
                    field(cells[3])
            });

    auto x = StaticMatrix<2, 4>({bp_.y, 1., 0., 0., bp_.x, 0., 1., 0.}) * A * b;

    return Vector2D(x(0, 0), x(1, 0));
}

Tensor2D GhostCellStencil::bpGrad(const VectorFiniteVolumeField &field) const
{
    auto cells = field.grid().findNearestNode(bp_).cells();

    Point2D x1 = cells[0].get().centroid();
    Point2D x2 = cells[1].get().centroid();
    Point2D x3 = cells[2].get().centroid();
    Point2D x4 = cells[3].get().centroid();

    auto A = inverse<4, 4>(
            {
                    x1.x * x1.y, x1.x, x1.y, 1.,
                    x2.x * x2.y, x2.x, x2.y, 1.,
                    x3.x * x3.y, x3.x, x3.y, 1.,
                    x4.x * x4.y, x4.x, x4.y, 1.,
            });

    auto b = StaticMatrix<4, 2>(
            {
                    field(cells[0]).x, field(cells[0]).y,
                    field(cells[1]).x, field(cells[1]).y,
                    field(cells[2]).x, field(cells[2]).y,
                    field(cells[3]).x, field(cells[3]).y
            });

    auto x = StaticMatrix<2, 4>({bp_.y, 1., 0., 0., bp_.x, 0., 1., 0.}) * A * b;

    return Tensor2D(x(0, 0), x(1, 0), x(1, 0), x(1, 1));
}