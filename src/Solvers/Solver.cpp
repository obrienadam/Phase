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

ScalarFiniteVolumeField& Solver::addScalarField(const Input& input, const std::string& name)
{
    typedef std::pair< std::string, ScalarFiniteVolumeField> Key;
    return (scalarFields_.insert(Key(name, ScalarFiniteVolumeField(input, grid_, name))).first)->second;
}

ScalarFiniteVolumeField& Solver::addScalarField(const std::string& name)
{
    typedef std::pair< std::string, ScalarFiniteVolumeField> Key;
    return (scalarFields_.insert(Key(name, ScalarFiniteVolumeField(grid_, name))).first)->second;
}

VectorFiniteVolumeField& Solver::addVectorField(const Input& input, const std::string& name)
{
    typedef std::pair< std::string, VectorFiniteVolumeField> Key;
    return (vectorFields_.insert(Key(name, VectorFiniteVolumeField(input, grid_, name))).first)->second;
}

VectorFiniteVolumeField& Solver::addVectorField(const std::string& name)
{
    typedef std::pair< std::string, VectorFiniteVolumeField> Key;
    return (vectorFields_.insert(Key(name, VectorFiniteVolumeField(grid_, name))).first)->second;
}
