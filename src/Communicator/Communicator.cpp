#include <cstdarg>

#include <mpi.h>

#include "Communicator.h"
#include "Exception.h"

MPI_Datatype Communicator::MPI_VECTOR2D_ = Communicator::initVector2DDataType();

Communicator::Communicator(MPI_Comm comm)
    :
      comm_(comm)
{
    scalarSendBuffers_.resize(nProcs());
    scalarRecvBuffers_.resize(nProcs());
    vector2DSendBuffers_.resize(nProcs());
    vector2DRecvBuffers_.resize(nProcs());

    recvIdOrdering_.resize(nProcs());
}

Communicator::~Communicator()
{

}

int Communicator::printf(const char *format, ...)
{
    int n = 0;

    if(mainProc())
    {
        va_list argsPtr;
        va_start(argsPtr, format);
        n = vfprintf(stdout, format, argsPtr);
        va_end(argsPtr);
    }

    return n;
}

//- Info
int Communicator::rank() const
{
    int rank;
    MPI_Comm_rank(comm_, &rank);

    return rank;
}

int Communicator::nProcs() const
{
    int nProcs;
    MPI_Comm_size(comm_, &nProcs);

    return nProcs;
}

//- Communication

void Communicator::broadcast(int root, std::vector<int> &ints) const
{
    MPI_Bcast(ints.data(), ints.size(), MPI_INT, root, comm_);
}

void Communicator::broadcast(int root, std::vector<double> &doubles) const
{
    MPI_Bcast(doubles.data(), doubles.size(), MPI_DOUBLE, root, comm_);
}

void Communicator::broadcast(int root, std::vector<Vector2D> &vector2Ds) const
{
    MPI_Bcast(vector2Ds.data(), vector2Ds.size(), MPI_VECTOR2D_, root, comm_);
}

void Communicator::scatter(int root, const std::vector<int> &send, std::vector<int> &recv)
{
    MPI_Scatter(send.data(), send.size()/nProcs(), MPI_INT, recv.data(), recv.size(), MPI_INT, root, comm_);
}

int Communicator::scatter(int root, const std::vector<int> &send) const
{
    int num[1];
    MPI_Scatter(send.data(), 1, MPI_INT, num, 1, MPI_INT, root, comm_);
    return num[0];
}

void Communicator::scatterv(int root, const std::vector<int>& sendBuff, const std::vector<int>& counts, std::vector<int>& recvBuff)
{
    std::vector<int> displ(nProcs(), 0);

    for(int i = 1; i < nProcs(); ++i)
        displ[i] = counts[i - 1];

    recvBuff.resize(counts[rank()]);
    MPI_Scatterv(sendBuff.data(), counts.data(), displ.data(), MPI_INT, recvBuff.data(), counts[rank()], MPI_INT, root, comm_);
}

//- Blocking point-to-point communication
void Communicator::send(int dest, const std::vector<double> &vals) const
{
    MPI_Send(vals.data(), vals.size(), MPI_DOUBLE, dest, MPI_ANY_TAG, comm_);
}

void Communicator::send(int dest, const std::vector<Vector2D> &vals) const
{
    MPI_Send(vals.data(), vals.size(), MPI_VECTOR2D_, dest, MPI_ANY_TAG, comm_);
}

void Communicator::recv(int source, std::vector<double> &vals) const
{
    MPI_Status status;
    MPI_Recv(vals.data(), vals.size(), MPI_DOUBLE, source, MPI_ANY_TAG, comm_, &status);
}

void Communicator::recv(int source, std::vector<Vector2D> &vals) const
{
    MPI_Status status;
    MPI_Recv(vals.data(), vals.size(), MPI_VECTOR2D_, source, MPI_ANY_TAG, comm_, &status);
}

//- Non-blocking point-to-point communication
void Communicator::isend(int dest, const::std::vector<double> &vals) const
{
    MPI_Request request;
    MPI_Isend(vals.data(), vals.size(), MPI_DOUBLE, dest, MPI_ANY_TAG, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::isend(int dest, const::std::vector<Vector2D>& vals) const
{
    MPI_Request request;
    MPI_Isend(vals.data(), vals.size(), MPI_VECTOR2D_, dest, MPI_ANY_TAG, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::ibsend(int dest, const std::vector<unsigned long> &vals) const
{
    MPI_Request request;
    MPI_Ibsend(vals.data(), vals.size(), MPI_UNSIGNED_LONG, dest, MPI_ANY_TAG, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::ibsend(int dest, const std::vector<double> &vals) const
{
    MPI_Request request;
    MPI_Ibsend(vals.data(), vals.size(), MPI_DOUBLE, dest, MPI_ANY_TAG, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::ibsend(int dest, const std::vector<Vector2D>& vals) const
{
    MPI_Request request;
    MPI_Ibsend(vals.data(), vals.size(), MPI_VECTOR2D_, dest, MPI_ANY_TAG, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<unsigned long> &vals) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_UNSIGNED_LONG, source, MPI_ANY_TAG, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<double> &vals) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_DOUBLE, source, MPI_ANY_TAG, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<Vector2D>& vals) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_VECTOR2D_, source, MPI_ANY_TAG, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::waitAll() const
{
    std::vector<MPI_Status> statuses(currentRequests_.size());
    MPI_Waitall(currentRequests_.size(), currentRequests_.data(), statuses.data());

    currentRequests_.clear();
}

//- Collective communication
double Communicator::sum(double val) const
{
    double result;
    MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return result;
}

double Communicator::min(double val) const
{
    double result;
    MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_MIN, comm_);
    return result;
}

double Communicator::max(double val) const
{
    double result;
    MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_MAX, comm_);
    return result;
}

//- Perform field communications across processes
void Communicator::sendMessages(ScalarFiniteVolumeField& field) const
{
    for(int procNo = 0; procNo < nProcs(); ++procNo)
    {
        scalarSendBuffers_[procNo].clear();

        for(const Cell& cell: sendBufferCells_[procNo])
            scalarSendBuffers_[procNo].push_back(field(cell));

        if(!scalarSendBuffers_[procNo].empty())
            isend(procNo, scalarSendBuffers_[procNo]);
    }

    for(int procNo = 0; procNo < nProcs(); ++procNo)
    {
        if(!scalarRecvBuffers_[procNo].empty())
            irecv(procNo, scalarRecvBuffers_[procNo]);
    }

    waitAll(); // Allow comms to finalize

    for(int procNo = 0; procNo < nProcs(); ++procNo) // Unload recv buffers
    {
        for(int i = 0; i < scalarRecvBuffers_.size(); ++i)
            field[recvIdOrdering_[procNo][i]] = scalarRecvBuffers_[procNo][i];
    }
}

void Communicator::sendMessages(VectorFiniteVolumeField& field) const
{
    for(int procNo = 0; procNo < nProcs(); ++procNo)
    {
        vector2DSendBuffers_[procNo].clear();

        for(const Cell& cell: sendBufferCells_[procNo])
            vector2DSendBuffers_[procNo].push_back(field(cell));

        if(!vector2DSendBuffers_[procNo].empty())
            isend(procNo, vector2DSendBuffers_[procNo]);
    }

    for(int procNo = 0; procNo < nProcs(); ++procNo)
    {
        if(!vector2DRecvBuffers_[procNo].empty())
            irecv(procNo, vector2DRecvBuffers_[procNo]);
    }

    waitAll(); // Allow comms to finalize

    for(int procNo = 0; procNo < nProcs(); ++procNo) // Unload recv buffers
    {
        for(int i = 0; i < vector2DRecvBuffers_.size(); ++i)
            field[recvIdOrdering_[procNo][i]] = vector2DRecvBuffers_[procNo][i];
    }
}

//- Private methods

MPI_Datatype Communicator::initVector2DDataType()
{
    MPI_Datatype MPI_VECTOR2D;

    MPI_Type_vector(1, 2, 2, MPI_DOUBLE, &MPI_VECTOR2D);
    MPI_Type_commit(&MPI_VECTOR2D);

    return MPI_VECTOR2D;
}
