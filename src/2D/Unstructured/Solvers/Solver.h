#ifndef PHASE_SOLVER_H
#define PHASE_SOLVER_H

#include "System/SolverInterface.h"

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"
#include "FiniteVolume/Field/TensorFiniteVolumeField.h"
#include "FiniteVolume/ImmersedBoundary/ImmersedBoundary.h"

class Solver : public SolverInterface
{
public:
    //- Constructors
    Solver(const Input &input);

    //- Info
    virtual std::string info() const
    { return ""; }

    //- Initialize
    virtual void initialize() {}

    //- Solve
    virtual Scalar solve(Scalar timeStep) = 0;

    //- Time step methods
    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const = 0;

    Scalar maxTimeStep() const
    { return maxTimeStep_; }

    int printf(const char *format, ...) const;

    const Communicator& comm() const
    { return grid_->comm(); }

    Scalar getStartTime(const Input &input) const;

    //- Field management

    std::shared_ptr<FiniteVolumeField<int>> addIntegerField(const std::string &name);

    std::shared_ptr<ScalarFiniteVolumeField> addScalarField(const Input &input, const std::string &name);

    std::shared_ptr<ScalarFiniteVolumeField> addScalarField(const std::string &name);

    std::shared_ptr<VectorFiniteVolumeField> addVectorField(const Input &input, const std::string &name);

    std::shared_ptr<VectorFiniteVolumeField> addVectorField(const std::string &name);

    std::shared_ptr<ScalarFiniteVolumeField> addScalarField(const std::shared_ptr<ScalarFiniteVolumeField>& field);

    std::shared_ptr<VectorFiniteVolumeField> addVectorField(const std::shared_ptr<VectorFiniteVolumeField>& field);

    std::shared_ptr<TensorFiniteVolumeField> addTensorField(const std::shared_ptr<TensorFiniteVolumeField>& field);

    //- Field data structures
    const std::unordered_map<std::string, std::shared_ptr<FiniteVolumeField<int>>> &integerFields() const
    { return integerFields_; }

    const std::unordered_map<std::string, std::shared_ptr<ScalarFiniteVolumeField>> &scalarFields() const
    { return scalarFields_; }

    const std::unordered_map<std::string, std::shared_ptr<VectorFiniteVolumeField>> &vectorFields() const
    { return vectorFields_; }

    //- Field lookup
    const std::shared_ptr<FiniteVolumeField<int>> &integerField(const std::string &name) const
    { return integerFields_.find(name)->second; }

    const std::shared_ptr<ScalarFiniteVolumeField> &scalarField(const std::string &name) const
    { return scalarFields_.find(name)->second; }

    const std::shared_ptr<VectorFiniteVolumeField> &vectorField(const std::string &name) const
    { return vectorFields_.find(name)->second; }

    //- Grid
    std::shared_ptr<FiniteVolumeGrid2D> grid()
    { return grid_; }

    std::shared_ptr<const FiniteVolumeGrid2D> grid() const
    { return grid_; }

    //- Immersed boundary
    std::shared_ptr<ImmersedBoundary> ib()
    { return ib_; }

    std::shared_ptr<const ImmersedBoundary> ib() const
    { return ib_; }

    //- ICs/IBs
    void setInitialConditions(const Input &input);

protected:

    void setCircle(const Circle &circle, Scalar innerValue, ScalarFiniteVolumeField &field);

    void setCircle(const Circle &circle, const Vector2D &innerValue, VectorFiniteVolumeField &field);

    void setCircleSector(const Circle &circle,
                         Scalar thetaMin,
                         Scalar thetaMax,
                         Scalar innerValue,
                         ScalarFiniteVolumeField &field);

    void setBox(const Polygon &box, Scalar innerValue, ScalarFiniteVolumeField &field);

    void setBox(const Polygon &box, const Vector2D &innerValue, VectorFiniteVolumeField &field);

    void setRotating(const std::string &function,
                     Scalar amplitude,
                     const Vector2D &center,
                     ScalarFiniteVolumeField &field);

    void setRotating(const std::string &xFunction, const std::string &yFunction, const Vector2D &amplitude,
                     const Vector2D &center, VectorFiniteVolumeField &field);

    virtual void restartSolution(const Input &input);

    std::shared_ptr<FiniteVolumeGrid2D> grid_;

    std::shared_ptr<ImmersedBoundary> ib_;

    //- Fields and geometries
    mutable std::unordered_map<std::string, std::shared_ptr<FiniteVolumeField<int>>> integerFields_;

    mutable std::unordered_map<std::string, std::shared_ptr<ScalarFiniteVolumeField>> scalarFields_;

    mutable std::unordered_map<std::string, std::shared_ptr<VectorFiniteVolumeField>> vectorFields_;

    mutable std::unordered_map<std::string, std::shared_ptr<TensorFiniteVolumeField>> tensorFields_;

    //- Solver parameters
    bool restartedSolution_;

    Scalar maxTimeStep_;
};

#endif
