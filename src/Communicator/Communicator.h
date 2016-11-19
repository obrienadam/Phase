#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

//-------------------------------------------------
// This class will be used for managing the future
// planned MPI support for Phase. It currently
// is not actually used in the code, but is placed
// here for development and testing.
//-------------------------------------------------

#include <mpi.h>

#include "FiniteVolumeGrid2D.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

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
    void broadcast(int root, std::string& str) const;
    void broadcast(int root, std::vector<Vector2D> &vector2Ds) const;

    //- Collectives
    int scatter(int root, const std::vector<int> &send) const;
    std::vector<unsigned long> allGather(unsigned long val) const;

    //- Blocking point-to-point communication
    void send(int dest, const std::vector<unsigned long>& vals, int tag = 0) const;
    void send(int dest, const std::vector<Vector2D>& vals, int tag = 0) const;

    void recv(int source, std::vector<double>& vals) const;
    void recv(int source, std::vector<Vector2D> &vals) const;

    //- Non-blocking point-to-point communication
    void ibsend(int dest, const std::vector<unsigned long>& vals, int tag = 0) const;

    void irecv(int source, std::vector<unsigned long>& vals, int tag = 0) const;
    void irecv(int source, std::vector<Vector2D>& vals, int tag = 0) const;

    void waitAll() const;

    //- Collective communications
    double sum(double val) const;
    double min(double val) const;
    double max(double val) const;

private:

    static MPI_Datatype MPI_VECTOR2D_;

    MPI_Comm comm_;

    Field<Label> globalCellIds_;

    mutable std::vector< std::vector<double> > scalarSendBuffers_, scalarRecvBuffers_;
    mutable std::vector< std::vector<Vector2D> > vector2DSendBuffers_, vector2DRecvBuffers_;

    std::vector<CellGroup> sendBufferCells_;
    std::vector<CellGroup> recvBufferCells_;

    std::vector< std::vector<Label> > recvIdOrdering_;

    mutable std::vector<MPI_Request> currentRequests_;
};

#endif
