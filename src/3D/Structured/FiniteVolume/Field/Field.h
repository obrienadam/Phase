#ifndef PHASE_FIELD_H
#define PHASE_FIELD_H

#include "Structured/StructuredGrid3D/StructuredGrid3D.h"

template<class T>
class Field
{
public:

    Field(const std::string &name,
          const std::shared_ptr<const StructuredGrid3D> &grid,
          bool isCellField = true,
          bool isFaceField = true);

    Field(const std::string &name,
          const Input &input,
          const std::shared_ptr<const StructuredGrid3D> &grid,
          bool isCellField = true,
          bool isFaceField = true);

    const std::string &name() const
    { return _name; }

    T& operator()(const Cell& cell)
    { return _cells[cell.id()]; }

    const T& operator()(const Cell& cell) const
    { return _cells[cell.id()]; }

    const std::vector<T> &cellData() const
    { return _cells; }

    //- Grid access
    const std::shared_ptr<const StructuredGrid3D> &grid() const
    { return _grid; }

protected:

    std::string _name;

    std::vector<T> _cells, _iFaces, _jFaces, _kFaces;

    std::shared_ptr<const StructuredGrid3D> _grid;
};

#include "Field.tpp"

#endif
