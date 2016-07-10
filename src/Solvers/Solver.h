#ifndef SOLVER_H
#define SOLVER_H

#include "FiniteVolumeGrid2D.h"
#include "Input.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "Circle.h"

class Solver
{
public:

    enum TimeDependent{ON, OFF};

    Solver(const FiniteVolumeGrid2D& grid, const Input& input);

    virtual std::string info();
    virtual Scalar solve(Scalar timeStep) = 0;
    virtual Scalar computeMaxTimeStep(Scalar maxCo) const = 0;

    ScalarFiniteVolumeField& addScalarField(const Input& input, const std::string& name);
    ScalarFiniteVolumeField& addScalarField(const std::string& name);

    VectorFiniteVolumeField& addVectorField(const Input& input, const std::string& name);
    VectorFiniteVolumeField& addVectorField(const std::string& name);

    std::vector<Polygon>& addGeometries(const std::string& name);

    const FiniteVolumeGrid2D& grid() const { return grid_; }
    std::map<std::string, ScalarFiniteVolumeField >& scalarFields() const { return scalarFields_; }
    std::map<std::string, VectorFiniteVolumeField >& vectorFields() const { return vectorFields_; }
    std::map<std::string, std::vector<Polygon> >& geometries() const { return geometries_; }

    void setInitialConditions(const Input& input);

protected:

    void setCircle(const Circle& circle, Scalar innerValue, ScalarFiniteVolumeField& field);
    void setCircle(const Circle& circle, const Vector2D& innerValue, VectorFiniteVolumeField& field);

    void setBox(const Polygon& box, Scalar innerValue, ScalarFiniteVolumeField& field);
    void setBox(const Polygon& box, const Vector2D& innerValue, VectorFiniteVolumeField& field);

    void setRotating(const std::string& function, Scalar amplitude, const Vector2D& center, ScalarFiniteVolumeField& field);
    void setRotating(const std::string& xFunction, const std::string& yFunction, const Vector2D& amplitude, const Vector2D& center, VectorFiniteVolumeField& field);

    const FiniteVolumeGrid2D& grid_;

    mutable std::map<std::string, ScalarFiniteVolumeField > scalarFields_;
    mutable std::map<std::string, VectorFiniteVolumeField > vectorFields_;
    mutable std::map<std::string, std::vector<Polygon> > geometries_;

    TimeDependent timeDependent_;
    int maxIterations_;
};

#endif
