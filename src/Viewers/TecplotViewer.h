#ifndef TECPLOT_VIEWER_H
#define TECPLOT_VIEWER_H

#include <mpi.h> // Necessary to get rid of missing def errors
#include <vector>
#include <fstream>

#include "Viewer.h"

class TecplotViewer : public Viewer
{
public:

    TecplotViewer(const Input& input, const Communicator& comm, const Solver& solver);
    ~TecplotViewer();

    void close();

    void write(Scalar solutionTime, const Communicator &comm);

private:
    std::vector<std::ofstream> outFiles_;
    int zoneNo_ = 1;
};

#endif
