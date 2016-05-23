#ifndef HEIGHT_FUNCTION_H
#define HEIGHT_FUNCTION_H

#include "SurfaceTensionForce.h"
#include "Polygon.h"
#include "ScalarFiniteVolumeField.h"

class HeightFunction : public SurfaceTensionForce
{
public:

    HeightFunction(const Input &input, const ScalarFiniteVolumeField& gamma);

    virtual VectorFiniteVolumeField compute();

private:

    Point2D loc_;
    std::vector<Polygon> cols_;
};

#endif
