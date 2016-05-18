#ifndef PLIC_H
#define PLIC_H

#include "Equation.h"

namespace plic
{

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, ScalarFiniteVolumeField &field);

Polygon computeInterfacePolygon(const Cell& cell, Scalar gamma, const Vector2D& normal);
}

#endif
