#ifndef PHASE_CGNS_VIEWER_H
#define PHASE_CGNS_VIEWER_H

#include <unordered_set>

#include "System/CgnsFile.h"
#include "System/Input.h"

#include "Structured/Solvers/Solver.h"

class CgnsViewer
{
public:

    CgnsViewer(const Input &input, const std::weak_ptr<const Solver> &solver);

    void write(Scalar timeStep);

protected:

    std::weak_ptr<const Solver> _solver;

    std::unordered_set<std::string> _scalarFieldNames;

};

#endif
