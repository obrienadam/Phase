#include <cstdarg>

#include <mpi.h>
#include <metis.h>

#include "Communicator.h"
#include "Exception.h"

MPI_Datatype Communicator::MPI_VECTOR2D_ = Communicator::initVector2DDataType();

Communicator::Communicator(MPI_Comm comm)
    :
      comm_(comm)
{

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

//- Partitioning
void Communicator::partitionGrid(const FiniteVolumeGrid2D &grid)
{
    using namespace std;

    int nElems = grid.cells().size();
    int nNodes = grid.nodes().size();

    vector<idx_t> elems;
    vector<idx_t> elemInds;

    elems.reserve(4*nElems);
    elemInds.reserve(nElems + 1);
    elemInds.push_back(0);

    for(const Cell &cell: grid.cells())
    {
        const Size nVerts = cell.nodes().size();
        elemInds.push_back(elemInds.back() + nVerts);

        for(const Node &node: cell.nodes())
            elems.push_back(node.id());
    }

    vector<idx_t> elemPart(nElems);
    vector<idx_t> nodePart(nNodes);

    idx_t nCommon = 1;
    idx_t nParts = nProcs();
    idx_t objVal;

    if(mainProc())
    {
        printf("Partitioning mesh into %d subdomains...\n", nParts);
        int status = METIS_PartMeshDual(&nElems, &nNodes,
                                        elemInds.data(), elems.data(),
                                        NULL, NULL,
                                        &nCommon, &nParts,
                                        NULL, NULL, &objVal,
                                        elemPart.data(), nodePart.data());

        if(status != METIS_OK)
            throw Exception("Communicator", "partitionGrid", "an error was encountered while partitioning the grid.");
    }

    broadcast(mainProcNo(), elemPart);
    broadcast(mainProcNo(), nodePart);
}

//- Private methods
MPI_Datatype Communicator::initVector2DDataType()
{
    MPI_Datatype MPI_VECTOR2D;

    MPI_Type_vector(1, 2, 2, MPI_DOUBLE, &MPI_VECTOR2D);
    MPI_Type_commit(&MPI_VECTOR2D);

    return MPI_VECTOR2D;
}
