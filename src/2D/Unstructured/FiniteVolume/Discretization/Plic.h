#ifndef PLIC_H
#define PLIC_H

#include "2D/Unstructured/FiniteVolume/Equation/FiniteVolumeEquation.h"

namespace plic {

FiniteVolumeEquation<Scalar> div(const VectorFiniteVolumeField &u,
                                 const VectorFiniteVolumeField &gradGamma,
                                 ScalarFiniteVolumeField &gamma,
                                 Scalar timeStep);

/*
FiniteVolumeEquation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField
&u, const VectorFiniteVolumeField& gradGamma, ScalarFiniteVolumeField &field,
Scalar timeStep, std::vector<Polygon>& plicPolygons);

Polygon computeInterfacePolygon(const Cell& cell, Scalar &gamma, const Vector2D&
normal); Polygon computeFluxPolygon(const BoundaryLink &link, const Vector2D&
uf, Scalar timeStep, int componentNo);
*/
} // namespace plic

#endif
