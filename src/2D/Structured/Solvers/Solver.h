#ifndef PHASE_SOLVER_H
#define PHASE_SOLVER_H

#include <unordered_set>
#include <functional>

#include "System/Input.h"
#include "System/SolverInterface.h"

#include "StructuredGrid2D/StructuredGrid2D.h"
#include "FiniteVolume/Field/ScalarField.h"
#include "FiniteVolume/Field/VectorField.h"

class Solver: public SolverInterface
{
public:

    Solver(const Input &input, const std::shared_ptr<StructuredGrid2D> &grid);

    virtual Scalar getStartTime() const override
    { return 0.; }

    virtual std::string info() const
    { }

    virtual void initialize() override
    { }

    virtual void setInitialConditions(const Input &input) override;

    virtual void setInitialConditions(const CommandLine &cl, const Input &input) override;

    //- Timestep
    virtual Scalar maxTimeStep() const
    { return _maxTimeStep; }

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar timeStep) const override
    { return _maxTimeStep; }

    virtual int printf(const char *format, ...) const override
    { }

    //- comm and grid access
    virtual const Communicator& comm() const override
    { return _grid->comm(); }

    std::shared_ptr<StructuredGrid2D> &grid()
    { return _grid; }

    std::shared_ptr<const StructuredGrid2D> grid() const
    { return _grid; }

    //- Field management

    template<class T>
    Field<T> &addField(const std::string &name);

protected:

    struct Hash
    {
        std::size_t operator ()(const std::shared_ptr<ScalarField> &field) const
        { return std::hash<std::string>{}(field->name()); }

        std::size_t operator ()(const std::shared_ptr<VectorField> &field) const
        { return std::hash<std::string>{}(field->name()); }
    };

    struct EqualTo
    {
        bool operator()(const std::shared_ptr<ScalarField> &lhs, const std::shared_ptr<ScalarField> &rhs) const
        { return lhs->name() == rhs->name(); }

        bool operator()(const std::shared_ptr<VectorField> &lhs, const std::shared_ptr<VectorField> &rhs) const
        { return lhs->name() == rhs->name(); }
    };

    //- Grid
    std::shared_ptr<StructuredGrid2D> _grid;

    //- Fields
    std::unordered_set<std::shared_ptr<ScalarField>, Hash, EqualTo> _scalarFields;

    std::unordered_set<std::shared_ptr<VectorField>, Hash, EqualTo> _vectorFields;

    //- Timestep
    Scalar _startTime, _maxTimeStep;
};


#endif
