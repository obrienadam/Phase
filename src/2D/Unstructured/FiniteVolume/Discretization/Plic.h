#ifndef PLIC_H
#define PLIC_H

#include "FiniteVolume/Equation/Equation.h"

namespace plic {

    Equation<Scalar> div(const VectorFiniteVolumeField& u,
                         const VectorFiniteVolumeField& gradGamma,
                         ScalarFiniteVolumeField& gamma,
                         Scalar timeStep);

/*
Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u, const VectorFiniteVolumeField& gradGamma, ScalarFiniteVolumeField &field, Scalar timeStep,
                                      std::vector<Polygon>& plicPolygons);

Polygon computeInterfacePolygon(const Cell& cell, Scalar &gamma, const Vector2D& normal);
Polygon computeFluxPolygon(const BoundaryLink &link, const Vector2D& uf, Scalar timeStep, int componentNo);
*/
}

#endif
