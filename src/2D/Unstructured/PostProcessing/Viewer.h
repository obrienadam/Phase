#ifndef PHASE_VIEWER_H
#define PHASE_VIEWER_H

#include "System/Input.h"
#include "System/Communicator.h"
#include "Solvers/Solver.h"

class Viewer
{
public:

    Viewer(const Input& input, const Solver& solver);

    virtual void write(Scalar solutionTime) = 0;

protected:

    const Solver& solver_;

    std::string filename_;

    std::unordered_set<std::string> integerFields_, scalarFields_, vectorFields_;

};

#include "CgnsViewer.h"

#endif
