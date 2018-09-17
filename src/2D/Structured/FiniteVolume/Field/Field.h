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

    //- History
    Field<T>& oldField(int i)
    { return *_prevFields[i].second; }

    const Field<T>& oldField(int i) const
    { return *_prevFields[i].second; }

    void savePreviousSolution(Scalar timeStep, int nSolutions);

    //- Parallel
    void sendMessages(bool sync = true);

protected:

    const std::string &_name;

    std::shared_ptr<const StructuredGrid2D> _grid;

    std::vector<T> _cellData, _faceData;

    //- Previous Fields
    std::vector<std::pair<Scalar, std::shared_ptr<Field<T>>>> _prevFields;

    //- Misc buffers for parallel comms
    std::vector<std::vector<T>> _sendBuffers;

    std::vector<std::vector<T>> _recvBuffers;
};

#include "Field.tpp"

#endif
