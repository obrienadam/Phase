#include "Field.h"

template<class T>
Field<T>::Field(const std::string &name,
                const std::shared_ptr<const StructuredGrid3D> &grid,
                bool isCellField,
                bool isFaceField)
    :
      _name(name),
      _grid(grid)
{
    if(isCellField)
        _cells.resize(_grid->nCells());

    if(isFaceField)
        _faces.resize(_grid->nFaces());
}

template<class T>
Field<T>::Field(const std::string &name,
                const Input &input,
                const std::shared_ptr<const StructuredGrid3D> &grid,
                bool isCellField,
                bool isFaceField)
    :
      Field(name, grid, isCellField, isFaceField)
{

}

template<class T>
const std::unique_ptr<BoundaryCondition<T>> &Field<T>::bc(const BoundaryPatch &patch) const
{
    auto it = _bcs.find(patch.name());
    return it != _bcs.end() ? it->second: nullptr;
}
