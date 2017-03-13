#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <mpi.h>
#include <vector>

#include "Vector2D.h"

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

    MPI_Comm communicator() const
    { return comm_; }

    //- Sync
    void barrier() const;

    //- Communication
    int broadcast(int root, int integer) const;

    void broadcast(int root, std::vector<int> &ints) const;

    void broadcast(int root, std::vector<unsigned long> &vals) const;

    void broadcast(int root, std::vector<double> &doubles) const;

    void broadcast(int root, std::vector<Vector2D> &vector2Ds) const;

    //- Collectives
    int scatter(int root, const std::vector<int> &send) const;

    std::vector<int> allGather(int val) const;

    std::vector<unsigned long> allGather(unsigned long val) const;

    //- Blocking point-to-point communication
    void ssend(int dest, const std::vector<int> &vals, int tag = 0) const;

    void ssend(int dest, const std::vector<unsigned long> &vals, int tag = 0) const;

    void ssend(int dest, const std::vector<double> &vals, int tag = 0) const;

    void ssend(int dest, const std::vector<Vector2D> &vals, int tag = 0) const;

    void recv(int source, std::vector<int> &vals) const;

    void recv(int source, std::vector<unsigned long> &vals) const;

    void recv(int source, std::vector<double> &vals) const;

    void recv(int source, std::vector<Vector2D> &vals) const;

    //- Non-blocking point-to-point communication
    void irecv(int source, std::vector<int> &vals, int tag = 0) const;

    void irecv(int source, std::vector<unsigned long> &vals, int tag = 0) const;

    void irecv(int source, std::vector<double> &vals, int tag = 0) const;

    void irecv(int source, std::vector<Vector2D> &vals, int tag = 0) const;

    void waitAll() const;

    //- Collective communications
    long sum(long val) const;

    unsigned long sum(unsigned long val) const;

    double sum(double val) const;

    int min(int val) const;

    double min(double val) const;

    double max(double val) const;

private:

    static MPI_Datatype MPI_VECTOR2D_;

    MPI_Comm comm_;
    mutable std::vector<MPI_Request> currentRequests_;
};

#endif
