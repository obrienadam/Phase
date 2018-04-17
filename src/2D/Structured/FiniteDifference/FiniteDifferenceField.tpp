#include "FiniteDifferenceField.h"

template<class T>
FiniteDifferenceField<T>::FiniteDifferenceField(const std::shared_ptr<const StructuredGrid2D> &grid)
        :
        _grid(grid),
        _nodes(grid->nNodes())
{

}