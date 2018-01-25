#ifndef QUADRATIC_IBM_STENCIL_H
#define QUADRATIC_IBM_STENCIL_H

#include <valarray>

#include "ImmersedBoundary.h"

class QuadraticIbmStencil
{
public:
    QuadraticIbmStencil(const InteriorLink &link,
                        const ImmersedBoundary &ib);

    const std::vector<Ref<const Cell>> &cells() const
    { return cells_; }

    const std::valarray<Scalar> &coeffs() const
    { return coeffs_; }

    const Vector2D &src() const
    { return src_; }

protected:

    void initQuadraticCoeffs(const Cell &stCell,
                             const Cell &cell,
                             const Cell &ibCell,
                             const ImmersedBoundaryObject &ibObj);

    void initQuadraticCoeffs(const Cell &stCell,
                             const Cell &cell,
                             const Cell &ibCell,
                             const ImmersedBoundaryObject &ibObjL,
                             const ImmersedBoundaryObject &ibObjR);

    void initLinearCoeffs(const Cell &stCell,
                          const Cell &ibCell,
                          const ImmersedBoundaryObject &ibObj);

    void initLinearCoeffs(const Cell &stCell,
                          const Cell &ibCell,
                          const ImmersedBoundaryObject &ibObjL,
                          const ImmersedBoundaryObject &ibObjR);

    std::vector<Ref<const Cell>> cells_;
    std::valarray<Scalar> coeffs_;
    Vector2D src_;
};


#endif
