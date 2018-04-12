#include "QuadraticInterpolator.h"

void QuadraticInterpolator::setPoint(const Point2D &pt)
{
    pt_ = pt;
    cells_ = grid_.lock()->globalCells().nearestItems(pt, 6);
    isValid_ = cells_.size() == 6;

    if (isValid_)
    {
        Point2D x[] = {
                cells_[0].get().centroid(),
                cells_[1].get().centroid(),
                cells_[2].get().centroid(),
                cells_[3].get().centroid(),
                cells_[4].get().centroid(),
                cells_[5].get().centroid()
        };

        A_ = inverse<6>({
                                x[0].x * x[0].x, x[0].y * x[0].y, x[0].x * x[0].y, x[0].x, x[0].y, 1.,
                                x[1].x * x[1].x, x[1].y * x[1].y, x[1].x * x[1].y, x[1].x, x[1].y, 1.,
                                x[2].x * x[2].x, x[2].y * x[2].y, x[2].x * x[2].y, x[2].x, x[2].y, 1.,
                                x[3].x * x[3].x, x[3].y * x[3].y, x[3].x * x[3].y, x[3].x, x[3].y, 1.,
                                x[4].x * x[4].x, x[4].y * x[4].y, x[4].x * x[4].y, x[4].x, x[4].y, 1.,
                                x[5].x * x[5].x, x[5].y * x[5].y, x[5].x * x[5].y, x[5].x, x[5].y, 1.,
                        });
    }
}

Scalar QuadraticInterpolator::operator()(const ScalarFiniteVolumeField &field) const
{
    auto x = StaticMatrix<1, 6>({pt_.x * pt_.x, pt_.y * pt_.y, pt_.x * pt_.y, pt_.x, pt_.y, 1.});
    auto b = StaticMatrix<6, 1>({
                                        field(cells_[0]),
                                        field(cells_[1]),
                                        field(cells_[2]),
                                        field(cells_[3]),
                                        field(cells_[4]),
                                        field(cells_[5]),
                                });

    return (x * A_ * b)(0, 0);
}

Vector2D QuadraticInterpolator::operator()(const VectorFiniteVolumeField &field) const
{
    auto x = StaticMatrix<1, 6>({pt_.x * pt_.x, pt_.y * pt_.y, pt_.x * pt_.y, pt_.x, pt_.y, 1.});

    auto b = StaticMatrix<6, 2>({
                                        field(cells_[0]).x, field(cells_[0]).y,
                                        field(cells_[1]).x, field(cells_[1]).y,
                                        field(cells_[2]).x, field(cells_[2]).y,
                                        field(cells_[3]).x, field(cells_[3]).y,
                                        field(cells_[4]).x, field(cells_[4]).y,
                                        field(cells_[5]).x, field(cells_[5]).y,
                                });

    auto u = x * A_ * b;

    return Vector2D(u(0, 0), u(0, 1));
}

Vector2D QuadraticInterpolator::grad(const ScalarFiniteVolumeField &field) const
{
    auto b = StaticMatrix<6, 1>(
            {
                    field(cells_[0]),
                    field(cells_[1]),
                    field(cells_[2]),
                    field(cells_[3]),
                    field(cells_[4]),
                    field(cells_[5])
            });

    auto x = StaticMatrix<2, 6>({2. * pt_.x, 0., pt_.y, 1., 0., 0., 0., 2. * pt_.y, pt_.x, 0., 1., 0.}) * A_ * b;

    return Vector2D(x(0, 0), x(1, 0));
}

Tensor2D QuadraticInterpolator::grad(const VectorFiniteVolumeField &field) const
{
    auto b = StaticMatrix<6, 2>(
            {
                    field(cells_[0]).x, field(cells_[0]).y,
                    field(cells_[1]).x, field(cells_[1]).y,
                    field(cells_[2]).x, field(cells_[2]).y,
                    field(cells_[3]).x, field(cells_[3]).y,
                    field(cells_[4]).x, field(cells_[4]).y,
                    field(cells_[5]).x, field(cells_[5]).y,
            });

    auto x = StaticMatrix<2, 6>({2. * pt_.x, 0., pt_.y, 1., 0., 0., 0., 2. * pt_.y, pt_.x, 0., 1., 0.}) * A_ * b;

    return Tensor2D(x(0, 0), x(1, 0), x(1, 0), x(1, 1));
}
