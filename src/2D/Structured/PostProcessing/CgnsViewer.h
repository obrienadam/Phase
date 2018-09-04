#ifndef PHASE_CGNS_VIEWER_H
#define PHASE_CGNS_VIEWER_H

#include <unordered_set>

#include "Structured/Solvers/Solver.h"

class CgnsViewer
{
public:

    CgnsViewer(const Input &input, const Solver &solver);

    void write(Scalar time);

protected:

    const Solver &_solver;

    std::unordered_set<std::string> _scalarFields, _vectorFields;

};

#endif
