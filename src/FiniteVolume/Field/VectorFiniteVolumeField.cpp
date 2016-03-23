#include "VectorFiniteVolumeField.h"
#include "FiniteVolumeGrid2D.h"

VectorFiniteVolumeField::VectorFiniteVolumeField(const FiniteVolumeGrid2D &grid, std::string name)
    :
      Field<Vector2D>::Field(grid.cells.size(), Vector2D(), name),
      faces_(grid.faces.size(), Vector2D())
{

}
