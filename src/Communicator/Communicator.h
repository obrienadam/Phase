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

    static void init() { MPI_Init(NULL, NULL); }
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

    //- Communication
    void broadcast(int root, std::vector<int>& ints) const;
    void broadcast(int root, std::vector<double>& doubles) const;
    void broadcast(int root, std::vector<Vector2D> &vector2Ds) const;

    //- Blocking point-to-point communication
    void send(int dest, const std::vector<double>& vals) const;
    void send(int dest, const std::vector<Vector2D> &vals) const;

    void recv(int source, std::vector<double>& vals) const;
    void recv(int source, std::vector<Vector2D> &vals) const;

    //- Non-blocking point-to-point communication
    void isend(int dest, const::std::vector<double>& vals) const;
    void isend(int dest, const::std::vector<Vector2D>& vals) const;

    void ibsend(int dest, const std::vector<double>& vals) const; // buffered sends
    void ibsend(int dest, const std::vector<Vector2D>& vals) const;

    void irecv(int source, std::vector<double>& vals) const;
    void irecv(int source, std::vector<Vector2D>& vals) const;

    void waitAll() const;

    //- Collective communications
    double sum(double val) const;
    double min(double val) const;
    double max(double val) const;

    //- Partitions the grid on the communicator
    void partitionGrid(const FiniteVolumeGrid2D& grid);

    //- Perform field communications across processes
    void sendMessages(ScalarFiniteVolumeField& field) const;
    void sendMessages(VectorFiniteVolumeField& field) const;

private:

    static MPI_Datatype MPI_VECTOR2D_;
    static MPI_Datatype initVector2DDataType();

    MPI_Comm comm_;

    Field<Label> globalCellIds_;

    mutable std::vector< std::vector<double> > scalarSendBuffers_, scalarRecvBuffers_;
    mutable std::vector< std::vector<Vector2D> > vector2DSendBuffers_, vector2DRecvBuffers_;

    std::vector<UniqueCellGroup> sendCellGroups_;
    std::vector<CellGroup> recvCellGroups_;

    std::vector< std::vector<Label> > recvIdOrdering_;

    mutable std::vector<MPI_Request> currentRequests_;
};

#endif
