#ifndef SOLVER_FACTORY_H
#define SOLVER_FACTORY_H

#include "Solver.h"

class SolverFactory
{
public:
    enum SolverType
    {
        POISSON,
        FRACTIONAL_STEP,
        FRACTIONAL_STEP_EULER_LAGRANGE,
        FRACTIONAL_STEP_DIRECT_FORCING,
        FRACTIONAL_STEP_MULTIPHASE,
        FRACTIONAL_STEP_BOUSSINESQ
    };

    static std::shared_ptr<Solver> create(SolverType type,
                                          const Input &input);

    static std::shared_ptr<Solver> create(std::string type,
                                          const Input &input);

    static std::shared_ptr<Solver> create(const Input &input);
};


#endif //PHASE_SOLVERFACTORY_H
