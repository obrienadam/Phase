#include "SolverFactory.h"
#include "Poisson.h"
#include "FractionalStep.h"
#include "FractionalStepEulerLagrange.h"
#include "FractionalStepDirectForcing.h"
#include "FractionalStepMultiphase.h"
#include "FractionalStepBoussinesq.h"

std::shared_ptr<Solver> SolverFactory::create(SolverType type,
                                              const Input &input,
                                              const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
{
    switch (type)
    {
        case POISSON:
            return std::make_shared<Poisson>(input, grid);
        case FRACTIONAL_STEP:
            return std::make_shared<FractionalStep>(input, grid);
        case FRACTIONAL_STEP_EULER_LAGRANGE:
            return std::make_shared<FractionalStepEulerLagrange>(input, grid);
        case FRACTIONAL_STEP_DIRECT_FORCING:
            return std::make_shared<FractionalStepDirectForcing>(input, grid);
        case FRACTIONAL_STEP_MULTIPHASE:
            return std::make_shared<FractionalStepMultiphase>(input, grid);
        case FRACTIONAL_STEP_BOUSSINESQ:
            return std::make_shared<FractionalStepBoussinesq>(input, grid);
    }

    return nullptr;
}

std::shared_ptr<Solver> SolverFactory::create(std::string type,
                                              const Input &input,
                                              const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
{
    std::transform(type.begin(), type.end(), type.begin(), [](unsigned char c)
    {
        return std::tolower(c);
    });

    if (type == "poisson")
        return create(POISSON, input, grid);
    else if (type == "fractional step")
        return create(FRACTIONAL_STEP, input, grid);
    else if (type == "fractional step euler lagrange")
        return create(FRACTIONAL_STEP_EULER_LAGRANGE, input, grid);
    else if (type == "fractional step direct forcing")
        return create(FRACTIONAL_STEP_DIRECT_FORCING, input, grid);
    else if (type == "fractional step multiphase")
        return create(FRACTIONAL_STEP_MULTIPHASE, input, grid);
    else if (type == "fractional step boussinesq")
        return create(FRACTIONAL_STEP_BOUSSINESQ, input, grid);

    throw Exception("SolverFactory", "create", "solver \"" + type + "\" is not a valid solver type.");
}

std::shared_ptr<Solver> SolverFactory::create(const Input &input,
                                              const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
{
    return create(input.caseInput().get<std::string>("Solver.type"), input, grid);
}