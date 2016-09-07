#ifndef TOTAL_VARIATION_DIMINISHING
#define TOTAL_VARIATION_DIMINISHING

#include "Equation.h"
#include "TensorFiniteVolumeField.h"

// TVD schemes are currenlty not functioning properly on unstructured meshes!

namespace tvd
{

enum Limiter{SUPERBEE, MINMOD, OSHER, MUSCL, CHARM};

Equation<VectorFiniteVolumeField> div(const VectorFiniteVolumeField& u, VectorFiniteVolumeField& field, Limiter limiter = MUSCL);

}

#endif
