#include "FiniteVolumeField.h"

template <class T>
FiniteVolumeField<T>::FiniteVolumeField(const FiniteVolumeGrid2D &grid, const T &initialValue, const std::string &name)
    :
      grid_(grid),
      Field<T>::Field(grid.nCellsI(), grid.nCellsJ(), initialValue, name),
      facesI_(grid.nNodesI(), grid.nCellsJ(), initialValue, name + "_facesI"),
      facesJ_(grid.nCellsI(), grid.nNodesJ(), initialValue, name + "_facesJ")
{
    faceInterpolationMethod_ = VOLUME_WEIGHTED;
}

template <class T>
FiniteVolumeField<T>::FiniteVolumeField(const FiniteVolumeGrid2D& grid, const std::string &name)
    :
      FiniteVolumeField(grid, T(), name)
{

}
