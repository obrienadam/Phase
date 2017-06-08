#include "FiniteVolumeGrid2D.h"

template<typename T>
std::vector<Label> FiniteVolumeGrid2D::getCellIds(const T &cells)
{
    std::vector<Label> ids(cells.size());

    std::transform(cells.begin(), cells.end(), ids.begin(),
                   [](const Cell &cell) -> Label { return cell.id(); });

    return ids;
}

template<typename T>
void FiniteVolumeGrid2D::sendMessages(const Communicator &comm, std::vector<T> &data) const
{
    std::vector<std::vector<T>> recvBuffers(neighbouringProcs_.size());

    //- Post recvs first (non-blocking)
    for (int i = 0; i < neighbouringProcs_.size(); ++i)
    {
        recvBuffers[i].resize(bufferCellZones_[i]->size());
        comm.irecv(neighbouringProcs_[i], recvBuffers[i], comm.rank());
    }

    //- Send data (blocking sends)
    std::vector<T> sendBuffer;
    for (int i = 0; i < neighbouringProcs_.size(); ++i)
    {

        sendBuffer.resize(sendCellGroups_[i]->size());
        std::transform(sendCellGroups_[i]->begin(), sendCellGroups_[i]->end(),
                       sendBuffer.begin(), [&data](const Cell &cell) { return data[cell.id()]; });

        comm.ssend(neighbouringProcs_[i], sendBuffer, neighbouringProcs_[i]);
    }

    comm.waitAll();

    //- Unload recv buffers
    for (int i = 0; i < neighbouringProcs_.size(); ++i)
    {
        int j = 0;
        for (const Cell &cell: *bufferCellZones_[i])
            data[cell.id()] = recvBuffers[i][j++];
    }
}
