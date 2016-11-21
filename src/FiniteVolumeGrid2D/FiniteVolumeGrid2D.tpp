template<typename T>
void FiniteVolumeGrid2D::sendMessages(const Communicator &comm, std::vector<T> &data) const
{
    std::vector<std::vector<T>> recvBuffers(neighbouringProcs_.size());

    //- Post recvs first (non-blocking)
    for(int i = 0; i < neighbouringProcs_.size(); ++i)
    {
        recvBuffers[i].resize(procRecvOrder_[i].size());
        comm.irecv(neighbouringProcs_[i], recvBuffers[i]);
    }

    //- Send data (blocking sends)
    std::vector<T> sendBuffer;
    for(int i = 0; i < neighbouringProcs_.size(); ++i)
    {
        sendBuffer.clear();
        for(Label id: procSendOrder_[i])
            sendBuffer.push_back(data[id]);

        comm.send(neighbouringProcs_[i], sendBuffer);
    }

    //- Make sure all recv operations are completed
    comm.waitAll();

    //- Unload recv buffers
    for(int i = 0; i < neighbouringProcs_.size(); ++i)
        for(int j = 0; j < recvBuffers[i].size(); ++j)
            data[procRecvOrder_[i][j]] = recvBuffers[i][j];
}
