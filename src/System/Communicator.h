#ifndef PHASE_COMMUNICATOR_H
#define PHASE_COMMUNICATOR_H

#include <mpi.h>
#include <numeric>
#include <vector>

#include "2D/Geometry/Tensor2D.h"
#include "2D/Geometry/Vector2D.h"

#include "3D/Geometry/Point3D.h"
#include "3D/Geometry/Tensor3D.h"

class Communicator {
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

  int mainProcNo() const { return 0; }

  bool isMainProc() const { return rank() == mainProcNo(); }

  const MPI_Comm &communicator() const { return comm_; }

  //- Sync
  void barrier() const;

  //- Broadcast
  template <class T> T broadcast(int root, T val) const {
    MPI_Bcast(&val, sizeof(T), MPI_BYTE, root, comm_);
    return val;
  }

  template <class T> void broadcast(int root, std::vector<T> &vals) const {
    MPI_Bcast(vals.data(), sizeof(T) * vals.size(), MPI_BYTE, root, comm_);
  }

  //- gather
  template <class T> std::vector<T> gather(int root, const T &val) const {
    std::vector<T> result(nProcs());
    MPI_Gather(&val, sizeof(T), MPI_BYTE, result.data(), sizeof(T), MPI_BYTE,
               root, comm_);
    return result;
  }

  //- gatherv
  template <class T>
  std::vector<T> gatherv(int root, const std::vector<T> &vals) const {
    std::vector<int> sizes = gather(root, (int)(sizeof(T) * vals.size()));
    std::vector<T> result(std::accumulate(sizes.begin(), sizes.end(), 0) /
                          sizeof(T));
    std::vector<int> displs(sizes.size(), 0);
    std::partial_sum(sizes.begin(), sizes.end() - 1, displs.begin() + 1);

    MPI_Gatherv(vals.data(), sizeof(T) * vals.size(), MPI_BYTE, result.data(),
                sizes.data(), displs.data(), MPI_BYTE, root, comm_);

    return result;
  }

  //- Allgather
  template <class T> std::vector<T> allGather(const T &val) const {
    std::vector<T> result(nProcs());
    MPI_Allgather(&val, sizeof(T), MPI_BYTE, result.data(), sizeof(T), MPI_BYTE,
                  comm_);
    return result;
  }

  //- Allgatherv
  template <class T>
  std::vector<T> allGatherv(const std::vector<T> &vals) const {
    std::vector<int> sizes = allGather((int)(sizeof(T) * vals.size()));
    std::vector<T> result(std::accumulate(sizes.begin(), sizes.end(), 0) /
                          sizeof(T));
    std::vector<int> displs(sizes.size(), 0);
    std::partial_sum(sizes.begin(), sizes.end() - 1, displs.begin() + 1);

    MPI_Allgatherv(vals.data(), sizeof(T) * vals.size(), MPI_BYTE,
                   result.data(), sizes.data(), displs.data(), MPI_BYTE, comm_);

    return result;
  }

  //- Blocking point-to-point communication
  template <class T>
  void ssend(int dest, const std::vector<T> &vals,
             int tag = MPI_ANY_TAG) const {
    MPI_Ssend(vals.data(), sizeof(T) * vals.size(), MPI_BYTE, dest, tag, comm_);
  }

  template <class T>
  void recv(int source, std::vector<T> &vals, int tag = MPI_ANY_TAG) const {
    MPI_Status status;
    MPI_Recv(vals.data(), sizeof(T) * vals.size(), MPI_BYTE, source, tag, comm_,
             &status);
  }

  //- Non-blocking point-to-point communication

  template <class T>
  void isend(int dest, std::vector<T> &vals, int tag = MPI_ANY_TAG) const {
    MPI_Request request;
    MPI_Isend(vals.data(), sizeof(T) * vals.size(), MPI_BYTE, dest, tag, comm_,
              &request);
    currentRequests_.push_back(request);
  }

  template <class T>
  void irecv(int source, std::vector<T> &vals, int tag = MPI_ANY_TAG) const {
    MPI_Request request;
    MPI_Irecv(vals.data(), sizeof(T) * vals.size(), MPI_BYTE, source, tag,
              comm_, &request);
    currentRequests_.push_back(request);
  }

  void waitAll() const;

  //- Dynamic

  template <typename T> int probeSize(int source, int tag = MPI_ANY_TAG) const;

  //- Collective communications
  long sum(long val) const;

  unsigned long sum(unsigned long val) const;

  double sum(double val) const;

  Vector2D sum(const Vector2D &val) const;

  Tensor2D sum(const Tensor2D &val) const;

  Vector3D sum(const Vector3D &val) const;

  Tensor3D sum(const Tensor3D &val) const;

  int min(int val) const;

  double min(double val) const;

  double max(double val) const;

  //- Additional operators

private:
  static MPI_Datatype MPI_VECTOR2D_, MPI_TENSOR2D_;

  MPI_Comm comm_;

  mutable std::vector<MPI_Request> currentRequests_;

  mutable std::vector<MPI_Status> statuses_;
};

template <class T>
const Communicator &operator<<(const Communicator &comm, const T &val) {
  if (comm.isMainProc())
    std::cout << val;

  return comm;
}

#endif
