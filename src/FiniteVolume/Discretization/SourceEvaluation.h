#ifndef SOURCE_EVALUATION_H
#define SOURCE_EVALUATION_H

#include "VectorFiniteVolumeField.h"
#include "CellGroup.h"

namespace fv
{
    VectorFiniteVolumeField source(const CellGroup& cells, VectorFiniteVolumeField field);
    VectorFiniteVolumeField source(VectorFiniteVolumeField field);
}

#endif
