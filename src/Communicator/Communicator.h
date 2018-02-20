#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <mpi.h>
#include <vector>

#include "Vector2D.h"
#include "Tensor2D.h"

class Communicator
{
public:

    static void init(int argc, char *argv[]);

    static void finalize();

    Communicator(MPI_Comm comm = MPI_COMM_WORLD);

    ~Communicator();

    //- Printing

    int printf(const char *format, ...) const;

    //- Info
    int rank() const;

    int nProcs() const;

    int mainProcNo() const
    { return 0; }

    bool isMainProc() const
    { return rank() == mainProcNo(); }

    const MPI_Comm &communicator() const
    { return comm_; }

    //- Sync
    void barrier() const;

    //- Communication
    int broadcast(int root, int integer) const;

    Vector2D broadcast(int root, Vector2D vec) const;

    double broadcast(int root, double val) const;

    void broadcast(int root, std::vector<int> &ints) const;

    void broadcast(int root, std::vector<unsigned long> &vals) const;

    void broadcast(int root, std::vector<double> &doubles) const;

    void broadcast(int root, std::vector<Vector2D> &vector2Ds) const;

    //- Collectives

    //- Scatter
    int scatter(int root, const std::vector<int> &send) const;

    //- Allgather
    std::vector<int> allGather(int val) const;

    std::vector<unsigned long> allGather(unsigned long val) const;

    std::vector<Vector2D> allGather(const Vector2D &val) const;

    //- Allgatherv
    std::vector<Scalar> allGatherv(const std::vector<double> &vals) const;

    std::vector<Vector2D> allGatherv(const std::vector<Vector2D> &vals) const;

    //- gather
    std::vector<int> gather(int root, int val) const;

    std::vector<unsigned long> gather(int root, unsigned long val) const;

    //- gatherv
    std::vector<double> gatherv(int root, const std::vector<double> &vals) const;

    std::vector<Vector2D> gatherv(int root, const std::vector<Vector2D> &vals) const;

    //- Blocking point-to-point communication
    void ssend(int dest, const std::vector<int> &vals, int tag = MPI_ANY_TAG) const;

    void ssend(int dest, const std::vector<unsigned long> &vals, int tag = MPI_ANY_TAG) const;

    void ssend(int dest, const std::vector<double> &vals, int tag = MPI_ANY_TAG) const;

    void ssend(int dest, const std::vector<Vector2D> &vals, int tag = MPI_ANY_TAG) const;

    void ssend(int dest, const std::vector<Tensor2D> &vals, int tag = MPI_ANY_TAG) const;

    void ssend(int dest, unsigned long val, int tag = MPI_ANY_TAG) const;

    void recv(int source, std::vector<unsigned long> &vals, int tag = MPI_ANY_TAG) const;

    void recv(int source, std::vector<Vector2D> &vals, int tag = MPI_ANY_TAG) const;

    //- Non-blocking point-to-point communication

    void isend(int dest, const std::vector<unsigned long> &vals, int tag = MPI_ANY_TAG) const;

    void irecv(int source, std::vector<int> &vals, int tag = MPI_ANY_TAG) const;

    void irecv(int source, std::vector<unsigned long> &vals, int tag = MPI_ANY_TAG) const;

    void irecv(int source, std::vector<double> &vals, int tag = MPI_ANY_TAG) const;

    void irecv(int source, std::vector<Vector2D> &vals, int tag = MPI_ANY_TAG) const;

    void irecv(int source, std::vector<Tensor2D> &vals, int tag = MPI_ANY_TAG) const;

    void irecv(int source, unsigned long &val, int tag = MPI_ANY_TAG) const;

    void waitAll() const;

    //- Dynamic

    template<typename T>
    int probeSize(int source, int tag = MPI_ANY_TAG) const;

    //- Collective communications
    long sum(long val) const;

    unsigned long sum(unsigned long val) const;

    double sum(double val) const;

    Vector2D sum(const Vector2D &val) const;

    int min(int val) const;

    double min(double val) const;

    double max(double val) const;

private:

    static MPI_Datatype MPI_VECTOR2D_, MPI_TENSOR2D_;

    MPI_Comm comm_;
    mutable std::vector<MPI_Request> currentRequests_;
};

#endif
