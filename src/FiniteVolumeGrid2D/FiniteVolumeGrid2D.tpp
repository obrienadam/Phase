#include "FiniteVolumeGrid2D.h"

template<typename T>
void FiniteVolumeGrid2D::sendMessages(const Communicator &comm, std::vector<T> &data) const
{
    std::vector<std::vector<T>> recvBuffers(neighbouringProcs_.size());

    //- Post recvs first (non-blocking)
    for(int i = 0; i < neighbouringProcs_.size(); ++i)
    {
        recvBuffers[i].reserve(procRecvOrder_[i].size());
        recvBuffers[i].resize(procRecvOrder_[i].size());
        comm.irecv(neighbouringProcs_[i], recvBuffers[i], comm.rank());
    }

    //- Send data (blocking sends)
    std::vector<T> sendBuffer;
    for(int i = 0; i < neighbouringProcs_.size(); ++i)
    {
        sendBuffer.clear();
        sendBuffer.reserve(procSendOrder_[i].size());
        for(Label id: procSendOrder_[i])
            sendBuffer.push_back(data[id]);

        comm.ssend(neighbouringProcs_[i], sendBuffer, neighbouringProcs_[i]);
    }

    comm.waitAll();

    //- Unload recv buffers
    for(int i = 0; i < neighbouringProcs_.size(); ++i)
        for(int j = 0; j < recvBuffers[i].size(); ++j)
            data[procRecvOrder_[i][j]] = recvBuffers[i][j];
}
