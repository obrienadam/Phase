#include "SolverFactory.h"
#include "Poisson.h"
#include "FractionalStep.h"
#include "FractionalStepQuadraticIbm.h"
#include "FractionalStepEulerLagrange.h"
#include "FractionalStepDirectForcing.h"
#include "FractionalStepMultiphase.h"
#include "FractionalStepBoussinesq.h"

std::shared_ptr<Solver> SolverFactory::create(SolverType type,
                                              const Input &input)
{
    switch (type)
    {
        case POISSON:
            return std::make_shared<Poisson>(input);
        case FRACTIONAL_STEP:
            return std::make_shared<FractionalStep>(input);
        case FRACTIONAL_STEP_QUADRATIC_IBM:
            return std::make_shared<FractionalStepQuadraticIbm>(input);
        case FRACTIONAL_STEP_EULER_LAGRANGE:
            return std::make_shared<FractionalStepEulerLagrange>(input);
        case FRACTIONAL_STEP_DIRECT_FORCING:
            return std::make_shared<FractionalStepDirectForcing>(input);
        case FRACTIONAL_STEP_MULTIPHASE:
            return std::make_shared<FractionalStepMultiphase>(input);
        case FRACTIONAL_STEP_BOUSSINESQ:
            return std::make_shared<FractionalStepBoussinesq>(input);
    }

    return nullptr;
}

std::shared_ptr<Solver> SolverFactory::create(std::string type,
                                              const Input &input)
{
    std::transform(type.begin(), type.end(), type.begin(), [](unsigned char c)
    {
        return std::tolower(c);
    });

    if (type == "poisson")
        return create(POISSON, input);
    else if (type == "fractional step")
        return create(FRACTIONAL_STEP, input);
    else if (type == "fractional step quadratic ibm")
        return create(FRACTIONAL_STEP_QUADRATIC_IBM, input);
    else if (type == "fractional step euler lagrange")
        return create(FRACTIONAL_STEP_EULER_LAGRANGE, input);
    else if (type == "fractional step direct forcing")
        return create(FRACTIONAL_STEP_DIRECT_FORCING, input);
    else if (type == "fractional step multiphase")
        return create(FRACTIONAL_STEP_MULTIPHASE, input);
    else if (type == "fractional step boussinesq")
        return create(FRACTIONAL_STEP_BOUSSINESQ, input);

    throw Exception("SolverFactory", "create", "solver \"" + type + "\" is not a valid solver type.");
}

std::shared_ptr<Solver> SolverFactory::create(const Input &input)
{
    return create(input.caseInput().get<std::string>("Solver.type"), input);
}