#include "SourceEvaluation.h"

namespace fv
{
    VectorFiniteVolumeField source(const CellGroup& cells, VectorFiniteVolumeField field)
    {
        for(const Cell &cell: cells)
            field(cell) *= cell.volume();

        return field;
    }

    VectorFiniteVolumeField source(VectorFiniteVolumeField field)
    {
        for (const Cell &cell: field.grid.cellZone("fluid"))
            field(cell) *= cell.volume();

        return field;
    }
}
