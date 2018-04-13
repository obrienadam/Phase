#ifndef PHASE_SOLVER_INTERFACE_H
#define PHASE_SOLVER_INTERFACE_H

#include "Types/Types.h"

#include "Communicator.h"
#include "CommandLine.h"
#include "Input.h"

class SolverInterface
{
public:

    virtual Scalar getStartTime() const = 0;

    virtual std::string info() const
    {}

    virtual void initialize() = 0;

    virtual void setInitialConditions(const Input &input) = 0;

    virtual void setInitialConditions(const CommandLine &cl, const Input &input) = 0;

    virtual Scalar maxTimeStep() const = 0;

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar timeStep) const = 0;

    virtual Scalar solve(Scalar timeStep) = 0;

    virtual int printf(const char *format, ...) const = 0;

    virtual const Communicator& comm() const = 0;
};

#endif
