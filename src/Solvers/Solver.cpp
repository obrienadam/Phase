#include <boost/algorithm/string.hpp>

#include "Solver.h"

Solver::Solver(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      grid_(grid)
{
    std::string timeDependentOpt = input.caseInput().get<std::string>("Solver.timeDependent");
    boost::to_lower(timeDependentOpt);

    timeDependent_ = timeDependentOpt == "on" ? ON : OFF;
}

std::string Solver::info()
{
    return "Solver info:\n"
           "Time dependent: " + std::string((timeDependent_ == ON) ? "On" : "Off") + "\n"
           "Max Iterations: " + std::to_string(maxIterations_) + "\n";
}
