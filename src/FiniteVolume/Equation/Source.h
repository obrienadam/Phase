#ifndef SOURCE_H
#define SOURCE_H

#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

namespace fv
{
    namespace src
    {

        ScalarFiniteVolumeField div(const VectorFiniteVolumeField &field);

        VectorFiniteVolumeField laplacian(const ScalarFiniteVolumeField &gamma,
                                          const VectorFiniteVolumeField &field);

        VectorFiniteVolumeField div(const ScalarFiniteVolumeField &rho,
                                    const VectorFiniteVolumeField &u,
                                    const VectorFiniteVolumeField &field);
    }

    template <class T>
    FiniteVolumeField<T> source(const FiniteVolumeField<T> &field, const CellGroup& group)
    {
        VectorFiniteVolumeField srcField(field.gridPtr(), field.name(), T(), false, false);

        for(const Cell& cell: group)
            srcField(cell) = field(cell) * cell.volume();

        return srcField;
    }
}

#endif
