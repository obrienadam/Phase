#ifndef PHASE_SOLVER_H
#define PHASE_SOLVER_H

#include "System/SolverInterface.h"
#include "System/CommandLine.h"

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"
#include "FiniteVolume/Field/TensorFiniteVolumeField.h"
#include "FiniteVolume/ImmersedBoundary/ImmersedBoundary.h"

class CgnsFile;

class Solver : public SolverInterface
{
public:
    //- Constructors
    Solver(const Input &input, const std::shared_ptr<const FiniteVolumeGrid2D> &grid);

    //- Info
    virtual std::string info() const
    { return "Unknown solver type"; }

    //- Initialize
    virtual void initialize() {}

    //- Solve
    virtual Scalar solve(Scalar timeStep) = 0;

    //- Time step methods
    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const = 0;

    Scalar maxTimeStep() const
    { return maxTimeStep_; }

    int printf(const char *format, ...) const;

    const Communicator &comm() const
    { return grid_->comm(); }

    Scalar getStartTime() const;

    //- Field management

    template<class T>
    std::shared_ptr<FiniteVolumeField<T>> addField(const std::string &name, const std::shared_ptr<CellGroup> &cells = nullptr);

    template<class T>
    std::shared_ptr<FiniteVolumeField<T>> addField(const Input &input, const std::string& name, const std::shared_ptr<CellGroup> &cells = nullptr);

    template<class T>
    std::shared_ptr<FiniteVolumeField<T>> addField(const std::shared_ptr<FiniteVolumeField<T>> &field);

    //- Field data structures
    const std::unordered_map<std::string, std::shared_ptr<FiniteVolumeField<int>>> &integerFields() const
    { return integerFields_; }

    const std::unordered_map<std::string, std::shared_ptr<ScalarFiniteVolumeField>> &scalarFields() const
    { return scalarFields_; }

    const std::unordered_map<std::string, std::shared_ptr<VectorFiniteVolumeField>> &vectorFields() const
    { return vectorFields_; }

    //- Field lookup
    std::shared_ptr<FiniteVolumeField<int> > integerField(const std::string &name) const;

    std::shared_ptr<ScalarFiniteVolumeField> scalarField(const std::string &name) const;

    std::shared_ptr<VectorFiniteVolumeField> vectorField(const std::string &name) const;

    //- Grid
    const std::shared_ptr<const FiniteVolumeGrid2D> &grid() const
    { return grid_; }

    //- ICs
    void setInitialConditions(const Input &input);

    void setInitialConditions(const CommandLine &cl, const Input &input);

    //- IBs
    virtual std::shared_ptr<const ImmersedBoundary> ib() const
    { return nullptr; }

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

    virtual void restartFromCgnsViewer(const Input &input);

    virtual void restartFromCompactCgnsViewer(const Input &input);

    virtual int readLatestCgnsFlowSolution(const CgnsFile& file);

    std::shared_ptr<const FiniteVolumeGrid2D> grid_;

    std::shared_ptr<IndexMap> scalarIndexMap_, vectorIndexMap_;

    //- Fields and geometries
    mutable std::unordered_map<std::string, std::shared_ptr<FiniteVolumeField<int>>> integerFields_;

    mutable std::unordered_map<std::string, std::shared_ptr<ScalarFiniteVolumeField>> scalarFields_;

    mutable std::unordered_map<std::string, std::shared_ptr<VectorFiniteVolumeField>> vectorFields_;

    mutable std::unordered_map<std::string, std::shared_ptr<TensorFiniteVolumeField>> tensorFields_;

    //- Solver parameters
    Scalar startTime_, maxTimeStep_;

    //- Misc
    bool isRestart_;
};

#endif
