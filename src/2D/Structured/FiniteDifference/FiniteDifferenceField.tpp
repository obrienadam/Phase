#include "FiniteDifferenceField.h"

template<class T>
FiniteDifferenceField<T>::FiniteDifferenceField(const std::shared_ptr<const StructuredGrid2D> &grid,
                                                const std::string &name)
        :
        _grid(grid),
        _name(name),
        _nodes(grid->nNodes())
{

}