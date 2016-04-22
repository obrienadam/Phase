#ifndef HEIGHT_FUNCTION_H
#define HEIGHT_FUNCTION_H

#include "Polygon.h"
#include "ScalarFiniteVolumeField.h"

class HeightFunction : public Polygon
{
public:
    HeightFunction() {}
    HeightFunction(const Point2D& loc, const Vector2D& normal, Scalar colWidth, Scalar colHeight);

    Scalar computeHeight(const ScalarFiniteVolumeField& field);

    const Point2D& loc() const { return loc_; }
    const Vector2D& norm() const { return unitNormal_; }
    const Vector2D& tan() const { return unitTangent_; }

private:

    Point2D loc_;
    Vector2D unitNormal_, unitTangent_;
    Scalar colWidth_, colHeight_;
};

//- External functions
std::vector<HeightFunction> heightFunctionSet(const Point2D& centerLoc, const Vector2D& normal, Scalar colWidth, Scalar colHeight, int nHeightFunctions);
Scalar computeCurvature(const std::vector<HeightFunction>& hfs, ScalarFiniteVolumeField& field);

#endif
