#ifndef REYNOLDS_NUMBER_H
#define REYNOLDS_NUMBER_H

#include <fstream>

#include "PostProcessingObject.h"
#include "Cell.h"

class ReynoldsNumber: public PostProcessingObject
{
public:

    ReynoldsNumber(const Solver& solver, const Cell& cell);

    void compute(Scalar time);

private:

    const Cell& cell_;

    std::ofstream fout_;
};


#endif
