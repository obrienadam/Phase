#include <cstdarg>
#include <numeric>

#include <mpi.h>

#include "Communicator.h"
#include "Exception.h"

MPI_Datatype Communicator::MPI_VECTOR2D_;
MPI_Datatype Communicator::MPI_TENSOR2D_;

void Communicator::init(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Type_vector(1, 2, 2, MPI_DOUBLE, &MPI_VECTOR2D_);
    MPI_Type_vector(1, 4, 4, MPI_DOUBLE, &MPI_TENSOR2D_);
    MPI_Type_commit(&MPI_VECTOR2D_);
    MPI_Type_commit(&MPI_TENSOR2D_);
}

void Communicator::finalize()
{
    MPI_Finalize();
}

Communicator::Communicator(MPI_Comm comm)
        :
        comm_(comm)
{

}

Communicator::~Communicator()
{

}

int Communicator::printf(const char *format, ...) const
{
    int n = 0;

    if (isMainProc())
    {
        va_list argsPtr;
        va_start(argsPtr, format);
        n = vfprintf(stdout, format, argsPtr);
        va_end(argsPtr);
    }

    return n;
}

int Communicator::printf(const std::string& format, ...) const
{
    int n = 0;

    if (isMainProc())
    {
        va_list argsPtr;
        va_start(argsPtr, format);
        n = vfprintf(stdout, format.c_str(), argsPtr);
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

Vector2D Communicator::broadcast(int root, Vector2D vec) const
{
    MPI_Bcast(&vec, 1, MPI_VECTOR2D_, root, comm_);
    return vec;
}

double Communicator::broadcast(int root, double val) const
{
    MPI_Bcast(&val, 1, MPI_DOUBLE, root, comm_);
    return val;
}

//- Communication

void Communicator::broadcast(int root, std::vector<int> &ints) const
{
    MPI_Bcast(ints.data(), ints.size(), MPI_INT, root, comm_);
}

void Communicator::broadcast(int root, std::vector<unsigned long> &vals) const
{
    MPI_Bcast(vals.data(), vals.size(), MPI_UNSIGNED_LONG, root, comm_);
}

void Communicator::broadcast(int root, std::vector<double> &doubles) const
{
    MPI_Bcast(doubles.data(), doubles.size(), MPI_DOUBLE, root, comm_);
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

std::vector<int> Communicator::allGather(int val) const
{
    std::vector<int> result(nProcs());
    MPI_Allgather(&val, 1, MPI_INT, result.data(), 1, MPI_INT, comm_);
    return result;
}

std::vector<unsigned long> Communicator::allGather(unsigned long val) const
{
    std::vector<unsigned long> result(nProcs());
    MPI_Allgather(&val, 1, MPI_UNSIGNED_LONG, result.data(), 1, MPI_UNSIGNED_LONG, comm_);
    return result;
}

std::vector<Vector2D> Communicator::allGather(const Vector2D &val) const
{
    std::vector<Vector2D> result(nProcs());
    MPI_Allgather(&val, 1, MPI_VECTOR2D_, result.data(), 1, MPI_VECTOR2D_, comm_);
    return result;
}

std::vector<Scalar> Communicator::allGatherv(const std::vector<double>& vals) const
{
    std::vector<int> sizes = allGather((int)vals.size());
    std::vector<double> result(std::accumulate(sizes.begin(), sizes.end(), 0));
    std::vector<int> displs(1, 0);
    std::partial_sum(sizes.begin(), sizes.end() - 1, std::back_inserter(displs));
    MPI_Allgatherv(vals.data(), vals.size(), MPI_DOUBLE, result.data(), sizes.data(), displs.data(), MPI_DOUBLE, comm_);

    return result;
}

std::vector<Vector2D> Communicator::allGatherv(const std::vector<Vector2D>& vals) const
{
    std::vector<int> sizes = allGather((int)vals.size());
    std::vector<Vector2D> result(std::accumulate(sizes.begin(), sizes.end(), 0));
    std::vector<int> displs(1, 0);
    std::partial_sum(sizes.begin(), sizes.end() - 1, std::back_inserter(displs));
    MPI_Allgatherv(vals.data(), vals.size(), MPI_VECTOR2D_, result.data(), sizes.data(), displs.data(), MPI_VECTOR2D_, comm_);

    return result;
}

std::vector<int> Communicator::gather(int root, int val) const
{
    std::vector<int> result(nProcs());
    MPI_Gather(&val, 1, MPI_INT, result.data(), 1, MPI_INT, root, comm_);
    return result;
}

std::vector<unsigned long> Communicator::gather(int root, unsigned long val) const
{
    std::vector<unsigned long> result(nProcs());
    MPI_Gather(&val, 1, MPI_UNSIGNED_LONG, result.data(), 1, MPI_UNSIGNED_LONG, root, comm_);
    return result;
}

std::vector<double> Communicator::gatherv(int root, const std::vector<double> &vals) const
{
    std::vector<int> sizes = gather(root, (int)vals.size());
    std::vector<double> result(std::accumulate(sizes.begin(), sizes.end(), 0));
    std::vector<int> displs(sizes.size(), 0);
    std::partial_sum(sizes.begin(), sizes.end() - 1, displs.begin() + 1);

    MPI_Gatherv(vals.data(), vals.size(), MPI_DOUBLE, result.data(), sizes.data(), displs.data(), MPI_DOUBLE, root, comm_);
    return result;
}

std::vector<Vector2D> Communicator::gatherv(int root, const std::vector<Vector2D>& vals) const
{
    std::vector<int> sizes = gather(root, (int)vals.size());
    std::vector<Vector2D> result(std::accumulate(sizes.begin(), sizes.end(), 0));
    std::vector<int> displs(sizes.size(), 0);
    std::partial_sum(sizes.begin(), sizes.end() - 1, displs.begin() + 1);

    MPI_Gatherv(vals.data(), vals.size(), MPI_VECTOR2D_, result.data(), sizes.data(), displs.data(), MPI_VECTOR2D_, root, comm_);
    return result;
}

void Communicator::ssend(int dest, const std::vector<int> &vals, int tag) const
{
    MPI_Ssend(vals.data(), vals.size(), MPI_INT, dest, tag, comm_);
}

void Communicator::ssend(int dest, const std::vector<unsigned long> &vals, int tag) const
{
    MPI_Ssend(vals.data(), vals.size(), MPI_UNSIGNED_LONG, dest, tag, comm_);
}

void Communicator::ssend(int dest, const std::vector<double> &vals, int tag) const
{
    MPI_Ssend(vals.data(), vals.size(), MPI_DOUBLE, dest, tag, comm_);
}

void Communicator::ssend(int dest, const std::vector<Vector2D> &vals, int tag) const
{
    MPI_Ssend(vals.data(), vals.size(), MPI_VECTOR2D_, dest, tag, comm_);
}

void Communicator::ssend(int dest, const std::vector<Tensor2D>& vals, int tag) const
{
    MPI_Ssend(vals.data(), vals.size(), MPI_TENSOR2D_, dest, tag, comm_);
}

void Communicator::ssend(int dest, unsigned long val, int tag) const
{
    MPI_Ssend(&val, 1, MPI_UNSIGNED_LONG, dest, tag, comm_);
}

void Communicator::recv(int source, std::vector<unsigned long> &vals, int tag) const
{
    MPI_Status status;
    MPI_Recv(vals.data(), vals.size(), MPI_UNSIGNED_LONG, source, tag, comm_, &status);
}

void Communicator::recv(int source, std::vector<Vector2D> &vals, int tag) const
{
    MPI_Status status;
    MPI_Recv(vals.data(), vals.size(), MPI_VECTOR2D_, source, tag, comm_, &status);
}

void Communicator::isend(int dest, const std::vector<unsigned long> &vals, int tag) const
{
    MPI_Request request;
    MPI_Isend(vals.data(), vals.size(), MPI_UNSIGNED_LONG, dest, tag, comm_, &request);
    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<int> &vals, int tag) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_INT, source, tag, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<unsigned long> &vals, int tag) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_UNSIGNED_LONG, source, tag, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<double> &vals, int tag) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_DOUBLE, source, tag, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<Vector2D> &vals, int tag) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_VECTOR2D_, source, tag, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, std::vector<Tensor2D>& vals, int tag) const
{
    MPI_Request request;
    MPI_Irecv(vals.data(), vals.size(), MPI_TENSOR2D_, source, tag, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::irecv(int source, unsigned long &val, int tag) const
{
    MPI_Request request;
    MPI_Irecv(&val, 1, MPI_UNSIGNED_LONG, source, tag, comm_, &request);

    currentRequests_.push_back(request);
}

void Communicator::waitAll() const
{
    std::vector<MPI_Status> statuses(currentRequests_.size());
    MPI_Waitall(currentRequests_.size(), currentRequests_.data(), statuses.data());

    currentRequests_.clear();
}

template<>
int Communicator::probeSize<unsigned long>(int source, int tag) const
{
    MPI_Status status;
    int count;
    MPI_Probe(source, tag, comm_, &status);
    MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
    return count;
}

template <>
int Communicator::probeSize<double>(int source, int tag) const
{
    MPI_Status status;
    int count;
    MPI_Probe(source, tag, comm_, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    return count;
}

template <>
int Communicator::probeSize<Vector2D>(int source, int tag) const
{
    MPI_Status status;
    int count;
    MPI_Probe(source, tag, comm_, &status);
    MPI_Get_count(&status, MPI_VECTOR2D_, &count);
    return count;
}

long Communicator::sum(long val) const
{
    long result;
    MPI_Allreduce(&val, &result, 1, MPI_LONG, MPI_SUM, comm_);
    return result;
}

unsigned long Communicator::sum(unsigned long val) const
{
    unsigned long result;
    MPI_Allreduce(&val, &result, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm_);
    return result;
}

//- Collective communication

double Communicator::sum(double val) const
{
    double result;
    MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return result;
}

Vector2D Communicator::sum(const Vector2D& val) const
{
    std::vector<Vector2D> vals = allGather(val);
    return std::accumulate(vals.begin(), vals.end(), Vector2D(0., 0.));
}

int Communicator::min(int val) const
{
    int result;
    MPI_Allreduce(&val, &result, 1, MPI_INT, MPI_MIN, comm_);
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

