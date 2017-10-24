#ifndef QUADRATIC_IBM_STENCIL_H
#define QUADRATIC_IBM_STENCIL_H

#include "ImmersedBoundary.h"

class QuadraticIbmStencil
{
public:
    QuadraticIbmStencil(const Cell &cell,
                        const Cell &ibCell,
                        const ImmersedBoundary &ib,
                        Scalar flux = 1.);

    const std::vector<Ref<const Cell>> &cells() const
    { return cells_; }

    const std::vector<Scalar> &coeffs() const
    { return coeffs_; }

    const Vector2D &src() const
    { return src_; }

protected:

    void initQuadraticCoeffs(const Cell &stCell,
                             const Cell &cell,
                             const Cell &ibCell,
                             const ImmersedBoundaryObject &ibObj);

    void initQuadraticCoeffs(const ImmersedBoundaryObject &ibObjL,
                             const Cell &stCell,
                             const Cell &cell,
                             const Cell &ibCell,
                             const ImmersedBoundaryObject &ibObjR);

    void initLinearCoeffs(const Cell &stCell,
                          const Cell &ibCell,
                          const ImmersedBoundaryObject &ibObj);

    void initLinearCoeffs(const ImmersedBoundaryObject &ibObjL,
                          const Cell &stCell,
                          const Cell &ibCell,
                          const ImmersedBoundaryObject &ibObjR);

    std::vector<Ref<const Cell>> cells_;
    std::vector<Scalar> coeffs_;
    Vector2D src_;
};


#endif
