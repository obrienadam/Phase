#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

//-------------------------------------------------
// This class will be used for managing the future
// planned MPI support for Phase. It currently
// is not actually used in the code, but is placed
// here for development and testing.
//-------------------------------------------------

#include <mpi.h>
#include <vector>

#include "Vector2D.h"

class Communicator
{
public:

    static void init(int argc, char* argv[]);
    static void finalize() { MPI_Finalize(); }

    Communicator(MPI_Comm comm);
    ~Communicator();

    //- Printing
    int printf(const char* format, ...);

    //- Info
    int rank() const;
    int nProcs() const;

    int mainProcNo() const { return 0; }
    bool mainProc() const { return rank() == mainProcNo(); }

    //- Sync
    void barrier() const;

    //- Communication
    int broadcast(int root, int integer) const;

    void broadcast(int root, std::vector<int>& ints) const;
    void broadcast(int root, std::vector<double>& doubles) const;
    void broadcast(int root, std::vector<Vector2D> &vector2Ds) const;

    //- Collectives
    int scatter(int root, const std::vector<int> &send) const;
    std::vector<unsigned long> allGather(unsigned long val) const;

    //- Blocking point-to-point communication
    void send(int dest, const std::vector<int>& vals, int tag = 0) const;
    void send(int dest, const std::vector<unsigned long>& vals, int tag = 0) const;
    void send(int dest, const std::vector<double>& vals, int tag = 0) const;
    void send(int dest, const std::vector<Vector2D>& vals, int tag = 0) const;

    void recv(int source, std::vector<int>& vals) const;
    void recv(int source, std::vector<unsigned long> &vals) const;
    void recv(int source, std::vector<double>& vals) const;
    void recv(int source, std::vector<Vector2D> &vals) const;

    //- Non-blocking point-to-point communication
    void ibsend(int dest, const std::vector<unsigned long>& vals, int tag = 0) const;

    void irecv(int source, std::vector<unsigned long>& vals, int tag = 0) const;
    void irecv(int source, std::vector<double>& vals, int tag = 0) const;
    void irecv(int source, std::vector<Vector2D>& vals, int tag = 0) const;

    void waitAll() const;

    //- Collective communications
    double sum(double val) const;
    double min(double val) const;
    double max(double val) const;

private:

    static MPI_Datatype MPI_VECTOR2D_;

    MPI_Comm comm_;
    mutable std::vector<MPI_Request> currentRequests_;
};

#endif
