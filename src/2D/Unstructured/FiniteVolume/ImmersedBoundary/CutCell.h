#ifndef CUT_CELL_H
#define CUT_CELL_H

#include "FiniteVolumeGrid2D/Cell/Cell.h"
#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"
#include "CutCellLink.h"
#include "ImmersedBoundaryObject.h"

class CutCell
{
public:

    CutCell(const Cell &cell, const ImmersedBoundaryObject& ibObj);

    const Cell& cell() const { return cell_; }

    const Polygon& solid() const { return solid_; }

    const Polygon& fluid() const { return fluid_; }

    Scalar fluidVolume() const { return fluid_.area(); }

    Scalar solidVolume() const { return solid_.area(); }

    Scalar totalVolume() const { return cell_.volume(); }

    const Point2D& fluidCentroid() const { return fluid_.centroid(); }

    const Point2D& solidCentroid() const { return solid_.centroid(); }

    Scalar alpha() const { return fluidVolume() / cell_.volume(); }

    bool hasSolidFraction() const { return solid_.isValid(); }

    bool hasFluidFraction() const { return fluid_.isValid(); }

    bool isSmall() const;

    const Vector2D& solidFaceNorm() const { return bFaceNorm_; }

    const LineSegment2D& bFace() const { return bFace_; }

    const std::vector<CutCellLink>& neighbours() const { return cutCellLinks_; }

    const std::vector<BoundaryLink>& boundaries() const { return cell_.boundaries(); }

    const std::vector<Ref<const CutCellLink>> neighbours(const CellGroup& cellGroup) const;

    bool intersectsIbObj() const { return ibObj_; }

    const ImmersedBoundaryObject& ibObj() const { return *ibObj_; }

    //- Operators
    operator const Cell&() const
    { return cell_; }

private:

    const Cell &cell_;
    const ImmersedBoundaryObject* ibObj_ = nullptr;

    Vector2D bFaceNorm_;
    LineSegment2D bFace_;
    Polygon solid_, fluid_;

    std::vector<CutCellLink> cutCellLinks_;
};

#endif
