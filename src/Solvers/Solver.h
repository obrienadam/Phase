#ifndef SOLVER_H
#define SOLVER_H

#include "FiniteVolumeGrid2D.h"
#include "Input.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class Solver
{
public:

    enum TimeDependent{ON, OFF};

    Solver(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual std::string info();
    virtual Scalar solve(Scalar timeStep) = 0;

    ScalarFiniteVolumeField& addScalarField(const Input& input, const std::string& name);
    ScalarFiniteVolumeField& addScalarField(const std::string& name);

    VectorFiniteVolumeField& addVectorField(const Input& input, const std::string& name);
    VectorFiniteVolumeField& addVectorField(const std::string& name);

    const FiniteVolumeGrid2D& grid() const { return grid_; }
    std::map<std::string, ScalarFiniteVolumeField >& scalarFields() const { return scalarFields_; }
    std::map<std::string, VectorFiniteVolumeField >& vectorFields() const { return vectorFields_; }

protected:

    const FiniteVolumeGrid2D& grid_;

    mutable std::map<std::string, ScalarFiniteVolumeField > scalarFields_;
    mutable std::map<std::string, VectorFiniteVolumeField > vectorFields_;

    TimeDependent timeDependent_;
    int maxIterations_;
};

#endif
