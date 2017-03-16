#ifndef SOURCE_H
#define SOURCE_H

#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

namespace source
{

    ScalarFiniteVolumeField div(const VectorFiniteVolumeField &field);

    VectorFiniteVolumeField laplacian(const ScalarFiniteVolumeField& gamma, const VectorFiniteVolumeField &field);

    VectorFiniteVolumeField div(const ScalarFiniteVolumeField& rho, const VectorFiniteVolumeField& u, const VectorFiniteVolumeField& field);
}

#endif
