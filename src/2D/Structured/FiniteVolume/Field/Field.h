#ifndef PHASE_FIELD_H
#define PHASE_FIELD_H

#include "StructuredGrid2D/StructuredGrid2D.h"

template<class T>
class Field
{
public:

    Field(const std::string &name,
          const std::shared_ptr<const StructuredGrid2D> &grid,
          bool cellField = true,
          bool faceField = true);

    //- Access

    const std::string &name() const
    { return _name; }

    const std::shared_ptr<const StructuredGrid2D> &grid() const
    { return _grid; }

    T& operator()(const Cell &cell)
    { return _cellData[cell.lid()]; }

    const T& operator ()(const Cell &cell) const
    { return _cellData[cell.lid()]; }

    T& operator()(const Face &face)
    { return _faceData[face.lid()]; }

    const T& operator ()(const Face &face) const
    { return _faceData[face.lid()]; }

protected:

    const std::string &_name;

    std::shared_ptr<const StructuredGrid2D> _grid;

    std::vector<T> _cellData, _faceData;
};

#include "Field.tpp"

#endif
