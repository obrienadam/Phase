#ifndef PHASE_SOURCE_H
#define PHASE_SOURCE_H

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"

namespace src
{
    ScalarFiniteVolumeField div(const VectorFiniteVolumeField &field);

    ScalarFiniteVolumeField div(const VectorFiniteVolumeField &field, const CellGroup &cells);

    ScalarFiniteVolumeField laplacian(Scalar gamma,
                                      const ScalarFiniteVolumeField &phi);

    ScalarFiniteVolumeField laplacian(const ScalarFiniteVolumeField& gamma,
                                      const ScalarFiniteVolumeField& phi);

    template <class T>
    FiniteVolumeField<T> src(const FiniteVolumeField<T> &field)
    {
        VectorFiniteVolumeField srcField(field.grid(), field.name(), T(), false, false);

        for(const Cell& cell: field.cells())
            srcField(cell) = field(cell) * cell.volume();

        return srcField;
    }
}

#endif
