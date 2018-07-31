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
    {
        _iFaces.resize(_grid->ifaces().size());
        _jFaces.resize(_grid->jfaces().size());
        _kFaces.resize(_grid->kfaces().size());
    }
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
