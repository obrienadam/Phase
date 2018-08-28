#include "Field.h"

template<class T>
Field<T>::Field(const std::string &name,
                const std::shared_ptr<const StructuredGrid2D> &grid,
                bool cellField,
                bool faceField)
    :
      _name(name),
      _grid(grid)
{
    if(cellField)
        _cellData.resize(_grid->cells().size());

    if(faceField)
        _faceData.resize(_grid->faces().size());
}
