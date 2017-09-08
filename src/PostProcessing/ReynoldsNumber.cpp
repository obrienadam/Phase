#include "ReynoldsNumber.h"

ReynoldsNumber::ReynoldsNumber(const Solver &solver, const Cell &cell)
        :
        PostProcessingObject(solver),
        cell_(cell)
{
    fout_.open("Re_cell_" + std::to_string(cell_.id()) + ".dat");
    fout_ << "time\tRe\n";
    fout_.close();
}

void ReynoldsNumber::compute(Scalar time)
{
    fout_.open("Re_cell_" + std::to_string(cell_.id()) + ".dat", std::ofstream::app);
    fout_ << time << "\t" << 0. << "\n";
    fout_.close();
}