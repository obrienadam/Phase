#ifndef IB_VIEWER_H
#define IB_VIEWER_H

#include <mpi.h> // Necessary to get rid of missing def errors
#include <vector>
#include <fstream>

#include "Viewer.h"

class IbViewer : public Viewer
{
public:

    IbViewer(const Input& input, const Solver& solver);
    ~IbViewer();

    void close();

    void write(Scalar solutionTime);

private:

    void writeForces(Scalar solutionTime) const;

    const ImmersedBoundary &ib_;
    std::vector<std::ofstream> ibFiles_;
    std::vector<std::ofstream> ibForces_;
    std::ofstream cutCellFile_;
    int zoneNo_ = 1;
};

#endif
