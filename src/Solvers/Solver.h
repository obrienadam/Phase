#ifndef SOLVER_H
#define SOLVER_H

#include "FiniteVolumeGrid2D.h"
#include "Input.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "SparseMatrixSolver.h"
#include "Circle.h"
#include "ImmersedBoundary.h"
#include "VolumeIntegrator.h"
#include "ForceIntegrator.h"

class Solver
{
public:

    enum TimeDependent{ON, OFF};

    //- Constructors
    Solver(const Input& input, const Communicator& comm, FiniteVolumeGrid2D &grid);

    //- Info
    virtual std::string info() const;

    //- Solve
    virtual Scalar solve(Scalar timeStep) = 0;

    //- Time step methods
    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const = 0;
    Scalar maxTimeStep() const { return maxTimeStep_; }

    //- Field management
    FiniteVolumeField<int>& addIntegerField(const std::string& name);

    ScalarFiniteVolumeField& addScalarField(const Input& input, const std::string& name);
    ScalarFiniteVolumeField& addScalarField(const std::string& name);

    VectorFiniteVolumeField& addVectorField(const Input& input, const std::string& name);
    VectorFiniteVolumeField& addVectorField(const std::string& name);

    FiniteVolumeField<int>& getIntegerField(const std::string& name) { return integerFields_.find(name)->second; }

    std::vector<Polygon>& addGeometries(const std::string& name);

    std::map<std::string, FiniteVolumeField<int> >& integerFields() const { return integerFields_; }
    std::map<std::string, ScalarFiniteVolumeField >& scalarFields() const { return scalarFields_; }
    std::map<std::string, VectorFiniteVolumeField >& vectorFields() const { return vectorFields_; }
    std::map<std::string, std::vector<Polygon> >& geometries() const { return geometries_; }

    //- Grid
    FiniteVolumeGrid2D& grid() { return grid_; }
    const FiniteVolumeGrid2D& grid() const { return grid_; }

    //- Comm
    const Communicator& comm() const { return comm_; }

    //- ICs/IBs
    void setInitialConditions(const Input& input);

    std::vector<Ref<const ImmersedBoundaryObject>> ibObjs() const { return ibObjManager_.ibObjs(); }

    //- Integrators
    const std::vector<VolumeIntegrator>& volumeIntegrators() const { return volumeIntegrators_; }
    std::vector<VolumeIntegrator>& volumeIntegrators() { return volumeIntegrators_; }

    const std::vector<ForceIntegrator>& forceIntegrators() const { return forceIntegrators_; }
    std::vector<ForceIntegrator>& forceIntegrators() { return forceIntegrators_; }

protected:

    void setCircle(const Circle& circle, Scalar innerValue, ScalarFiniteVolumeField& field);
    void setCircle(const Circle& circle, const Vector2D& innerValue, VectorFiniteVolumeField& field);

    void setCircleSector(const Circle& circle, Scalar thetaMin, Scalar thetaMax, Scalar innerValue, ScalarFiniteVolumeField& field);

    void setBox(const Polygon& box, Scalar innerValue, ScalarFiniteVolumeField& field);
    void setBox(const Polygon& box, const Vector2D& innerValue, VectorFiniteVolumeField& field);

    void setRotating(const std::string& function, Scalar amplitude, const Vector2D& center, ScalarFiniteVolumeField& field);
    void setRotating(const std::string& xFunction, const std::string& yFunction, const Vector2D& amplitude, const Vector2D& center, VectorFiniteVolumeField& field);

    FiniteVolumeGrid2D& grid_;
    const Communicator comm_;

    //- Fields and geometries
    mutable std::map<std::string, FiniteVolumeField<int> > integerFields_;
    mutable std::map<std::string, ScalarFiniteVolumeField > scalarFields_;
    mutable std::map<std::string, VectorFiniteVolumeField > vectorFields_;
    mutable std::map<std::string, std::vector<Polygon> > geometries_;

    TimeDependent timeDependent_;
    Scalar timeStepRelaxation_, maxTimeStep_;

    ImmersedBoundary ibObjManager_;

    std::vector<VolumeIntegrator> volumeIntegrators_;
    std::vector<ForceIntegrator> forceIntegrators_;
};

#endif
