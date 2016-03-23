#include <string>

#include "ConstructSolver.h"
#include "Poisson.h"
#include "Simple.h"
#include "Exception.h"

std::shared_ptr<Solver> constructSolver(const FiniteVolumeGrid2D &grid, const Input &input)
{
    using namespace std;

    string solverType = input.caseInput().get<string>("Solver.type");

    if (solverType == "Poisson")
    {
        return shared_ptr<Poisson>(new Poisson(grid, input));
    }
    else if (solverType == "Simple")
    {
        return shared_ptr<Simple>(new Simple(grid, input));
    }
    else
    {
        throw Exception("", "constructSolver", "invalid solver type \"" + solverType + "\".");
    }
}
