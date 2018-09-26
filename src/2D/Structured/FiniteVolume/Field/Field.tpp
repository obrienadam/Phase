#include "System/Exception.h"

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

    _recvBuffers.resize(_grid->comm().nProcs());
    _sendBuffers.resize(_grid->comm().nProcs());
}

template<class T>
Field<T>::Field(const std::string &name,
                const std::shared_ptr<const StructuredGrid2D> &grid,
                const Input &input,
                bool cellField,
                bool faceField)
    :
      Field(name, grid, cellField, faceField)
{
    auto bcinput = input.boundaryInput().get_child_optional("Boundaries." + _name);

    if(!bcinput)
        throw Exception("Field<T>", "Field", "Missing boundary info for field \"" + _name + "\".");

    std::unordered_map<std::string, Coordinates::Direction> pmap = {
        {"x+", Coordinates::I_POS},
        {"x-", Coordinates::I_NEG},
        {"y+", Coordinates::J_POS},
        {"y-", Coordinates::J_NEG}
    };

    std::unordered_map<std::string, BoundaryCondition::Type> bcmap = {
        {"fixed", BoundaryCondition::FIXED},
        {"normal_gradient", BoundaryCondition::ZERO_GRADIENT}
    };


    for(const std::string &pname: {"x+", "x-", "y+", "y-"})
    {
        auto type = bcinput.get().get<std::string>(pname);
        _bctypes.emplace(pmap.at(pname), bcmap.at(type));
    }
}

template<class T>
void Field<T>::savePreviousSolution(Scalar timeStep, int nSolutions)
{
    auto soln = std::make_shared<Field<T>>(*this);
    soln->_prevFields.clear();
    _prevFields.insert(_prevFields.begin(), std::make_pair(timeStep, soln));
    _prevFields.resize(nSolutions);
}

template<class T>
void Field<T>::sendMessages(bool sync)
{
    //- Post non-blocking comms
    for(int proc = 0; proc < _grid->comm().nProcs(); ++proc)
    {
        //- Post the recv
        auto &recvBuffer = _recvBuffers[proc];
        recvBuffer.resize(_grid->recvBuffer(proc).size());

        _grid->comm().irecv(proc, recvBuffer, proc);

        //- Post the send
        auto &sendBuffer = _sendBuffers[proc];
        sendBuffer.resize(_grid->sendBuffer(proc).size());

        std::transform(_grid->sendBuffer(proc).begin(),
                       _grid->sendBuffer(proc).end(),
                       sendBuffer.begin(),
                       [this](const Cell &cell) { return _cellData[cell.lid()]; });

        _grid->comm().isend(proc, sendBuffer, _grid->comm().rank());
    }

    if(sync)
        _grid->comm().waitAll();
}
