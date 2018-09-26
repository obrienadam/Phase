#include <cstdarg>
#include <numeric>

#include <mpi.h>

#include "Communicator.h"

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

//- Communication

void Communicator::waitAll() const
{
    statuses_.resize(currentRequests_.size());
    MPI_Waitall(currentRequests_.size(), currentRequests_.data(), statuses_.data());
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

template<>
int Communicator::probeSize<double>(int source, int tag) const
{
    MPI_Status status;
    int count;
    MPI_Probe(source, tag, comm_, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &count);
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

Vector2D Communicator::sum(const Vector2D &val) const
{
    Vector2D result;
    MPI_Allreduce(&val, &result, 2, MPI_DOUBLE, MPI_SUM, comm_);
    return result;
}

Tensor2D Communicator::sum(const Tensor2D &val) const
{
    Tensor2D result;
    MPI_Allreduce(&val, &result, 4, MPI_DOUBLE, MPI_SUM, comm_);
    return result;
}

Vector3D Communicator::sum(const Vector3D &val) const
{
    Vector3D result;
    MPI_Allreduce(&val, &result, 3, MPI_DOUBLE, MPI_SUM, comm_);
    return result;
}

Tensor3D Communicator::sum(const Tensor3D &val) const
{
    Tensor3D result;
    MPI_Allreduce(&val, &result, 9, MPI_DOUBLE, MPI_SUM, comm_);
    return result;
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
