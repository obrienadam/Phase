#ifndef CICSAM_H
#define CICSAM_H

#include "Equation.h"
#include "ImmersedBoundaryObject.h"

namespace cicsam
{

enum Type{HC, UQ};

Scalar hc(Scalar gammaTilde, Scalar coD);
Scalar uq(Scalar gammaTilde, Scalar coD);

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u,
                                      const VectorFiniteVolumeField& gradField,
                                      ScalarFiniteVolumeField &field,
                                      Scalar timeStep,
                                      Type type = HC);

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u,
                                      const VectorFiniteVolumeField& gradField,
                                      const std::vector<ImmersedBoundaryObject> &ibObj,
                                      ScalarFiniteVolumeField &field,
                                      Scalar timeStep,
                                      Type type = HC);

Equation<ScalarFiniteVolumeField> div(const VectorFiniteVolumeField &u,
                                      const VectorFiniteVolumeField& gradField,
                                      const VectorFiniteVolumeField &m,
                                      ScalarFiniteVolumeField &field,
                                      Scalar timeStep);

}

#endif
