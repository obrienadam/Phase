#include "FiniteVolumeGrid2D.h"

template<class T>
void FiniteVolumeGrid2D::sendMessages(std::vector<T> &data) const
{
    if(!comm_ || comm_->nProcs() == 1)
        return;

    std::vector<std::vector<T>> recvBuffers(comm_->nProcs());

    //- Post recvs first (non-blocking)
    for(int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        if(bufferCellZones_[proc].empty())
            continue;

        recvBuffers[proc].resize(bufferCellZones_[proc].size());
        comm_->irecv(proc, recvBuffers[proc], proc);
    }

    //- Send data (blocking sends)
    std::vector<T> sendBuffer;
    for(int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        if(sendCellGroups_[proc].empty())
            continue;

        sendBuffer.resize(sendCellGroups_[proc].size());

        std::transform(sendCellGroups_[proc].begin(),
                       sendCellGroups_[proc].end(),
                       sendBuffer.begin(),
                       [&data](const Cell& cell) { return data[cell.id()]; });

        comm_->ssend(proc, sendBuffer, comm_->rank());
    }
    comm_->waitAll();

    //- Unload recv buffers
    for(int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        if(recvBuffers[proc].empty())
            continue;

        int i = 0;
        for(const Cell& cell: bufferCellZones_[proc])
            data[cell.id()] = recvBuffers[proc][i++];
    }
}

template<class T>
void FiniteVolumeGrid2D::sendMessages(std::vector<T> &data, Size nSets) const
{
    if(!comm_ || comm_->nProcs() == 1)
        return;

    std::vector<std::vector<T>> recvBuffers(comm_->nProcs());

    //- Post recvs first (non-blocking)
    for(int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        if(bufferCellZones_[proc].empty())
            continue;

        recvBuffers[proc].resize(bufferCellZones_[proc].size() * nSets);
        comm_->irecv(proc, recvBuffers[proc], proc);
    }

    //- Send data (blocking sends)
    std::vector<T> sendBuffer;
    for(int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        if(sendCellGroups_[proc].empty())
            continue;

        sendBuffer.resize(sendCellGroups_[proc].size() * nSets);

        Size i = 0;
        for(Size set = 0; set < nSets; ++set)
            for(const Cell& cell: sendCellGroups_[proc])
                sendBuffer[i++] = data[cell.id() + set * nCells()];

        comm_->ssend(proc, sendBuffer, comm_->rank());
    }
    comm_->waitAll();

    //- Unload recv buffers
    for(int proc = 0; proc < comm_->nProcs(); ++proc)
    {
        if(recvBuffers[proc].empty())
            continue;

        int i = 0;
        for(Size set = 0; set < nSets; ++set)
            for(const Cell& cell: bufferCellZones_[proc])
                data[cell.id() + set * nCells()] = recvBuffers[proc][i++];
    }
}
