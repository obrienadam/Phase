#ifndef QUADRATIC_INTERPOLATOR_H
#define QUADRATIC_INTERPOLATOR_H

#include "FiniteVolumeGrid2D.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "StaticMatrix.h"

class QuadraticInterpolator
{
public:

    QuadraticInterpolator(const std::weak_ptr<const FiniteVolumeGrid2D>& grid) : grid_(grid) {}

    void setPoint(const Point2D& pt);

    Scalar operator()(const ScalarFiniteVolumeField &field) const;

    Vector2D operator()(const VectorFiniteVolumeField &field) const;

    Vector2D grad(const ScalarFiniteVolumeField &field) const;

    Tensor2D grad(const VectorFiniteVolumeField &field) const;

    bool isValid() const
    { return isValid_; }

    const Point2D &point() const
    { return pt_; }

private:

    bool isValid_ = false;

    Point2D pt_;

    StaticMatrix<6, 6> A_;

    std::vector<Ref<const Cell>> cells_;

    std::weak_ptr<const FiniteVolumeGrid2D> grid_;
};


#endif
