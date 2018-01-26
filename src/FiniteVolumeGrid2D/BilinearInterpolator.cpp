#include "BilinearInterpolator.h"

void BilinearInterpolator::setPoint(const Point2D &pt)
{
    pt_ = pt;
    cells_ = grid_.lock()->findNearestNode(pt_).cells();
    isValid_ = cells_.size() == 4;

    if(isValid_)
    {
        Point2D x1 = cells_[0].get().centroid();
        Point2D x2 = cells_[1].get().centroid();
        Point2D x3 = cells_[2].get().centroid();
        Point2D x4 = cells_[3].get().centroid();

        A_ = inverse<4, 4>({
                                   x1.x * x1.y, x1.x, x1.y, 1.,
                                   x2.x * x2.y, x2.x, x2.y, 1.,
                                   x3.x * x3.y, x3.x, x3.y, 1.,
                                   x4.x * x4.y, x4.x, x4.y, 1.,
                           });
    }
}

Scalar BilinearInterpolator::operator()(const ScalarFiniteVolumeField &field) const
{
    auto x = StaticMatrix<1, 4>({pt_.x * pt_.y, pt_.x, pt_.y, 1.});
    auto b = StaticMatrix<4, 1>({
                                        field(cells_[0]),
                                        field(cells_[1]),
                                        field(cells_[2]),
                                        field(cells_[3])
                                });

    return (x * A_ * b)(0, 0);
}

Vector2D BilinearInterpolator::operator()(const VectorFiniteVolumeField &field) const
{
    auto x = StaticMatrix<1, 4>({pt_.x * pt_.y, pt_.x, pt_.y, 1.});

    auto b = StaticMatrix<4, 2>({
                                        field(cells_[0]).x, field(cells_[0]).y,
                                        field(cells_[1]).x, field(cells_[1]).y,
                                        field(cells_[2]).x, field(cells_[2]).y,
                                        field(cells_[3]).x, field(cells_[3]).y
                                });

    StaticMatrix<1, 2> u = x * A_ * b;

    return Vector2D(u(0, 0), u(0, 1));
}

Vector2D BilinearInterpolator::grad(const ScalarFiniteVolumeField &field) const
{
    auto b = StaticMatrix<4, 1>(
            {
                    field(cells_[0]),
                    field(cells_[1]),
                    field(cells_[2]),
                    field(cells_[3])
            });

    auto x = StaticMatrix<2, 4>({pt_.y, 1., 0., 0., pt_.x, 0., 1., 0.}) * A_ * b;

    return Vector2D(x(0, 0), x(1, 0));
}

Tensor2D BilinearInterpolator::grad(const VectorFiniteVolumeField &field) const
{
    auto b = StaticMatrix<4, 2>(
            {
                    field(cells_[0]).x, field(cells_[0]).y,
                    field(cells_[1]).x, field(cells_[1]).y,
                    field(cells_[2]).x, field(cells_[2]).y,
                    field(cells_[3]).x, field(cells_[3]).y
            });

    auto x = StaticMatrix<2, 4>({pt_.y, 1., 0., 0., pt_.x, 0., 1., 0.}) * A_ * b;

    return Tensor2D(x(0, 0), x(1, 0), x(1, 0), x(1, 1));
}