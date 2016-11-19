#include <cstdarg>

#include <mpi.h>

#include "Communicator.h"
#include "Exception.h"

MPI_Datatype Communicator::MPI_VECTOR2D_;

void Communicator::init(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Type_vector(1, 2, 2, MPI_DOUBLE, &MPI_VECTOR2D_);
    MPI_Type_commit(&MPI_VECTOR2D_);
}

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

void Communicator::barrier() const
{
    MPI_Barrier(comm_);
}

int Communicator::broadcast(int root, int integer) const
{
    MPI_Bcast(&integer, 1, MPI_INT, root, comm_);
    return integer;
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

void Communicator::broadcast(int root, std::string &str) const
{
    int size = broadcast(root, str.size());
    char cstr[size];
    strcpy(cstr, str.c_str());

    MPI_Bcast(cstr, size, MPI_CHAR, root, comm_);
    str = std::string(cstr);
}

void Communicator::broadcast(int root, std::vector<Vector2D> &vector2Ds) const
{
    MPI_Bcast(vector2Ds.data(), vector2Ds.size(), MPI_VECTOR2D_, root, comm_);
}

int Communicator::scatter(int root, const std::vector<int> &send) const
{
    int num[1];
    MPI_Scatter(send.data(), 1, MPI_INT, num, 1, MPI_INT, root, comm_);
    return num[0];
}

std::vector<unsigned long> Communicator::allGather(unsigned long val) const
{
    std::vector<unsigned long> result(nProcs());
    MPI_Allgather(&val, 1, MPI_UNSIGNED_LONG, result.data(), 1, MPI_UNSIGNED_LONG, comm_);

    return result;
}

void Communicator::send(int dest, const std::vector<unsigned long> &vals, int tag) const
{
    MPI_Send(vals.data(), vals.size(), MPI_UNSIGNED_LONG, dest, tag, comm_);
}

void Communicator::send(int dest, const std::vector<Vector2D> &vals, int tag) const
{
    MPI_Send(vals.data(), vals.size(), MPI_VECTOR2D_, dest, tag, comm_);
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
void Communicator::ibsend(int dest, const std::vector<unsigned long> &vals, int tag) const
{
    MPI_Request request;
    MPI_Ibsend(vals.data(), vals.size(), MPI_UNSIGNED_LONG, dest, tag, comm_, &request);
    //MPI_Bs

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<unsigned long> &vals, int tag) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_UNSIGNED_LONG, source, tag, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<Vector2D>& vals, int tag) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_VECTOR2D_, source, tag, comm_, &request);

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

