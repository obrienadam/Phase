#include "SourceEvaluation.h"

#include "FaceInterpolation.h"

namespace fv
{
    VectorFiniteVolumeField source(VectorFiniteVolumeField field)
    {
        for (const Cell &cell: field.grid.cellZone("fluid"))
            field(cell) *= cell.volume();

        return field;
    }
}
