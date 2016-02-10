#include "Solver.h"

Solver::Solver(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      grid_(grid)
{
    timeDependent_ = input.get<std::string>("Solver.timeDependent") == "ON" ? ON : OFF;
    maxIterations_ = input.get<int>("Solver.maxIterations");
}

std::string Solver::info()
{
    return "Solver info:\n"
           "Time dependent: " + std::string((timeDependent_ == ON) ? "On" : "Off") + "\n"
           "Max Iterations: " + std::to_string(maxIterations_) + "\n";
}
