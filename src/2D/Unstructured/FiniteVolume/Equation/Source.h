#ifndef SOURCE_H
#define SOURCE_H

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

    VectorFiniteVolumeField ftc(const ScalarFiniteVolumeField& cellWeight,
                                const ScalarFiniteVolumeField& faceWeight,
                                const VectorFiniteVolumeField& field,
                                const CellGroup& cells);

    VectorFiniteVolumeField ftc(const ScalarFiniteVolumeField& weight,
                                const VectorFiniteVolumeField& field,
                                const CellGroup& cells);

    template <class T>
    FiniteVolumeField<T> src(const FiniteVolumeField<T> &field, const CellGroup& group)
    {
        VectorFiniteVolumeField srcField(field.grid(), field.name(), T(), false, false);

        for(const Cell& cell: group)
            srcField(cell) = field(cell) * cell.volume();

        return srcField;
    }
}

#endif
