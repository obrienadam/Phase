#include <algorithm>

#include "CollisionModel.h"
#include "FiniteVolumeGrid2D/FiniteVolumeGrid2D.h"
#include "FiniteVolumeGrid2D/Face/FaceGroup.h"

CollisionModel::CollisionModel(Scalar eps, Scalar range)
{
    eps_ = eps;
    range_ = range;
}

Vector2D CollisionModel::force(const ImmersedBoundaryObject &ibObjP, const ImmersedBoundaryObject &ibObjQ) const
{
    if (ibObjP.shape().type() == Shape2D::CIRCLE && ibObjQ.shape().type() == Shape2D::CIRCLE)
    {
        const Circle &c1 = static_cast<const Circle &>(ibObjP.shape());
        const Circle &c2 = static_cast<const Circle &>(ibObjQ.shape());

        const Vector2D &xp = c1.centroid();
        const Vector2D &xq = c2.centroid();

        Scalar r1 = c1.radius();
        Scalar r2 = c2.radius();

        Scalar d = (xp - xq).mag();

        return d > r1 + r2 + range_ ? Vector2D(0., 0.) : (xp - xq) / eps_ * std::pow(r1 + r2 + range_ - d, 2);
    }
    else
        throw Exception("CollisionModel", "force", "unsupported shape type.");
}

Vector2D CollisionModel::force(const ImmersedBoundaryObject &ibObj, const FiniteVolumeGrid2D &grid) const
{
    Vector2D fc = Vector2D(0., 0.);

    if (ibObj.shape().type() == Shape2D::CIRCLE)
    {
        const Circle &c = static_cast<const Circle &>(ibObj.shape());
        const Vector2D &xp = c.centroid();
        Scalar r = c.radius();

        for (const FaceGroup &p: grid.patches())
            for (const Face &f: p.itemsCoveredBy(Circle(c.centroid(), r + range_)))
            {
                if(!grid.localCells().isInSet(f.lCell()))
                    continue;

                const Vector2D &xq = f.centroid();
                Scalar d = (xp - xq).mag();

                fc += (xp - xq) / eps_ * pow(r + range_ - d, 2);
            }
    }

    return fc;
}
