#include "GhostCellStencil.h"

GhostCellStencil::GhostCellStencil(const Cell &cell,
                                   const ImmersedBoundaryObject &ibObj)

        :
        cell_(std::cref(cell))
{
    bp_ = ibObj.nearestIntersect(cell.centroid());
    ip_ = 2. * bp_ - cell.centroid();
    nw_ = ibObj.nearestEdgeNormal(bp_).unitVec();
    dirichletCells_ = ibObj.grid()->findNearestNode(ip_).cells();
    neumannCells_ = ibObj.grid()->findNearestNode(bp_).cells();

    if (dirichletCells_.size() != 4 || neumannCells_.size() != 4)
    {
        std::ostringstream sout;
        sout << "number of image point cells must be 4. Boundary point = " << bp_ << ", image point = " << ip_ << ".";
        throw Exception("GhostCellStencil", "GhostCellStencil", sout.str());
    }

    initDirichletCoeffs();
    initNeumannCoeffs();
}

GhostCellStencil::GhostCellStencil(const Cell &cell,
                                   const Point2D &bp,
                                   const Vector2D &nw,
                                   const FiniteVolumeGrid2D &grid)
        :
        cell_(std::cref(cell))
{
    bp_ = bp;
    ip_ = 2 * bp - cell.centroid();
    nw_ = nw.unitVec();
    dirichletCells_ = grid.findNearestNode(ip_).cells();
    neumannCells_ = grid.findNearestNode(bp_).cells();

    if (dirichletCells_.size() != 4 || neumannCells_.size() != 4)
    {
        std::ostringstream sout;
        sout << "number of image point cells must be 4. Boundary point = " << bp_ << ", image point = " << ip_ << ".";
        throw Exception("GhostCellStencil", "GhostCellStencil", sout.str());
    }

    initDirichletCoeffs();
    initNeumannCoeffs();
}

Scalar GhostCellStencil::ipValue(const ScalarFiniteVolumeField &field) const
{
    auto c = Ad_ * StaticMatrix<4, 1>(
            {
                    field(dirichletCells_[0]),
                    field(dirichletCells_[1]),
                    field(dirichletCells_[2]),
                    field(dirichletCells_[3])
            });

    return c(0, 0) * ip_.x * ip_.y + c(1, 0) * ip_.x + c(2, 0) * ip_.y + c(3, 0);
}

Vector2D GhostCellStencil::ipValue(const VectorFiniteVolumeField &field) const
{
    auto c = Ad_ * StaticMatrix<4, 2>(
            {
                    field(dirichletCells_[0]).x, field(dirichletCells_[0]).y,
                    field(dirichletCells_[1]).x, field(dirichletCells_[1]).y,
                    field(dirichletCells_[2]).x, field(dirichletCells_[2]).y,
                    field(dirichletCells_[3]).x, field(dirichletCells_[3]).y
            });

    return Vector2D(
            c(0, 0) * ip_.x * ip_.y + c(1, 0) * ip_.x + c(2, 0) * ip_.y + c(3, 0),
            c(0, 1) * ip_.x * ip_.y + c(1, 1) * ip_.x + c(2, 1) * ip_.y + c(3, 1)
    );
}

Scalar GhostCellStencil::bpValue(const ScalarFiniteVolumeField &field) const
{
    auto c = An_ * StaticMatrix<4, 1>(
            {
                    field(neumannCells_[0]),
                    field(neumannCells_[1]),
                    field(neumannCells_[2]),
                    field(neumannCells_[3])
            });

    return c(0, 0) * bp_.x * bp_.y + c(1, 0) * bp_.x + c(2, 0) * bp_.y + c(3, 0);
}

Vector2D GhostCellStencil::bpValue(const VectorFiniteVolumeField &field) const
{
    auto c = An_ * StaticMatrix<4, 2>(
            {
                    field(neumannCells_[0]).x, field(neumannCells_[0]).y,
                    field(neumannCells_[1]).x, field(neumannCells_[1]).y,
                    field(neumannCells_[2]).x, field(neumannCells_[2]).y,
                    field(neumannCells_[3]).x, field(neumannCells_[3]).y
            });

    return Vector2D(
            c(0, 0) * bp_.x * bp_.y + c(1, 0) * bp_.x + c(2, 0) * bp_.y + c(3, 0),
            c(0, 1) * bp_.x * bp_.y + c(1, 1) * bp_.x + c(2, 1) * bp_.y + c(3, 1)
    );
}

Vector2D GhostCellStencil::bpGrad(const ScalarFiniteVolumeField &field) const
{
    auto b = StaticMatrix<4, 1>(
            {
                    field(neumannCells_[0]),
                    field(neumannCells_[1]),
                    field(neumannCells_[2]),
                    field(neumannCells_[3])
            });

    auto x = StaticMatrix<2, 4>({bp_.y, 1., 0., 0., bp_.x, 0., 1., 0.}) * An_ * b;

    return Vector2D(x(0, 0), x(1, 0));
}

Tensor2D GhostCellStencil::bpGrad(const VectorFiniteVolumeField &field) const
{
    auto b = StaticMatrix<4, 2>(
            {
                    field(neumannCells_[0]).x, field(neumannCells_[0]).y,
                    field(neumannCells_[1]).x, field(neumannCells_[1]).y,
                    field(neumannCells_[2]).x, field(neumannCells_[2]).y,
                    field(neumannCells_[3]).x, field(neumannCells_[3]).y
            });

    auto x = StaticMatrix<2, 4>({bp_.y, 1., 0., 0., bp_.x, 0., 1., 0.}) * An_ * b;
    return Tensor2D(x(0, 0), x(1, 0), x(1, 0), x(1, 1));
}

//- Private helpers

void GhostCellStencil::initDirichletCoeffs()
{
    Point2D x1 = dirichletCells_[0].get().centroid();
    Point2D x2 = dirichletCells_[1].get().centroid();
    Point2D x3 = dirichletCells_[2].get().centroid();
    Point2D x4 = dirichletCells_[3].get().centroid();

    bool ghostCellInStencil = false;
    for (const Cell &cell: dirichletCells_)
        if (cell_.get().id() == cell.id())
        {
            ghostCellInStencil = true;
            break;
        }

    Ad_ = inverse<4, 4>(
            {
                    x1.x * x1.y, x1.x, x1.y, 1.,
                    x2.x * x2.y, x2.x, x2.y, 1.,
                    x3.x * x3.y, x3.x, x3.y, 1.,
                    x4.x * x4.y, x4.x, x4.y, 1.,
            });

    if (ghostCellInStencil)
    {
        auto xd = StaticMatrix<1, 4>({bp_.x * bp_.y, bp_.x, bp_.y, 1.}) * Ad_;
        dirichletCoeffs_.assign(xd.data(), xd.data() + 4);
    }
    else
    {
        dirichletCells_.push_back(cell_);
        auto xd = StaticMatrix<1, 4>({ip_.x * ip_.y, ip_.x, ip_.y, 1.}) * Ad_ / 2.;
        dirichletCoeffs_.assign(xd.data(), xd.data() + 4);
        dirichletCoeffs_.push_back(1. / 2.);
    }
}

void GhostCellStencil::initNeumannCoeffs()
{
    Point2D x1 = neumannCells_[0].get().centroid();
    Point2D x2 = neumannCells_[1].get().centroid();
    Point2D x3 = neumannCells_[2].get().centroid();
    Point2D x4 = neumannCells_[3].get().centroid();

    bool ghostCellInStencil = false;
    for (const Cell &cell: neumannCells_)
        if (cell_.get().id() == cell.id())
        {
            ghostCellInStencil = true;
            break;
        }

    An_ = inverse<4, 4>(
            {
                    x1.x * x1.y, x1.x, x1.y, 1.,
                    x2.x * x2.y, x2.x, x2.y, 1.,
                    x3.x * x3.y, x3.x, x3.y, 1.,
                    x4.x * x4.y, x4.x, x4.y, 1.,
            });

    if (ghostCellInStencil)
    {
        auto xn = StaticMatrix<1, 4>({bp_.y * nw_.x + bp_.x * nw_.y, nw_.x, nw_.y, 0.}) * An_;
        neumannCoeffs_.assign(xn.data(), xn.data() + 4);
    }
    else
    {
        neumannCells_.push_back(cell_);
        auto xn = StaticMatrix<1, 4>({bp_.x * bp_.y, bp_.x, bp_.y, 1.}) * An_ / -length();
        neumannCoeffs_.assign(xn.data(), xn.data() + 4);
        neumannCoeffs_.push_back(1. / length());
    }
}