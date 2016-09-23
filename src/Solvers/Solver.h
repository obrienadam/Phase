#ifndef SOLVER_H
#define SOLVER_H

#include "FiniteVolumeGrid2D.h"
#include "Input.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "Circle.h"
#include "ImmersedBoundary.h"
#include "VolumeIntegrator.h"
#include "ForceIntegrator.h"

class Solver
{
public:

    enum TimeDependent{ON, OFF};

    Solver(const Input& input, FiniteVolumeGrid2D &grid);

    virtual std::string info() const;
    virtual Scalar solve(Scalar timeStep) = 0;

    virtual Scalar computeMaxTimeStep(Scalar maxCo, Scalar prevTimeStep) const = 0;
    Scalar maxTimeStep() const { return maxTimeStep_; }

    ScalarFiniteVolumeField& addScalarField(const Input& input, const std::string& name);
    ScalarFiniteVolumeField& addScalarField(const std::string& name);

    VectorFiniteVolumeField& addVectorField(const Input& input, const std::string& name);
    VectorFiniteVolumeField& addVectorField(const std::string& name);

    std::vector<Polygon>& addGeometries(const std::string& name);

    FiniteVolumeGrid2D& grid() { return grid_; }
    const FiniteVolumeGrid2D& grid() const { return grid_; }

    std::map<std::string, ScalarFiniteVolumeField >& scalarFields() const { return scalarFields_; }
    std::map<std::string, VectorFiniteVolumeField >& vectorFields() const { return vectorFields_; }
    std::map<std::string, std::vector<Polygon> >& geometries() const { return geometries_; }

    void setInitialConditions(const Input& input);

    const ImmersedBoundary& ib() const { return ib_; }

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

    mutable std::map<std::string, ScalarFiniteVolumeField > scalarFields_;
    mutable std::map<std::string, VectorFiniteVolumeField > vectorFields_;
    mutable std::map<std::string, std::vector<Polygon> > geometries_;

    TimeDependent timeDependent_;
    Scalar timeStepRelaxation_, maxTimeStep_;

    ImmersedBoundary ib_;

    std::vector<VolumeIntegrator> volumeIntegrators_;
    std::vector<ForceIntegrator> forceIntegrators_;
};

#endif
