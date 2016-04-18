#ifndef INPUT_H
#define INPUT_H

#include <string>

#include <boost/property_tree/ptree.hpp>

#include "Circle.h"

class Solver;
class ScalarFiniteVolumeField;
class VectorFiniteVolumeField;

class Input
{
public:

    Input(const std::string& caseDirectory = "case", const std::string& outputPath = "solution");

    void parseInputFile();

    std::string caseDirectory, outputPath;

    const boost::property_tree::ptree& caseInput() const { return caseInput_; }
    const boost::property_tree::ptree& boundaryInput() const { return boundaryInput_; }
    const boost::property_tree::ptree& initialConditionInput() const { return initialConditionInput_; }

    void setInitialConditions(const Solver &solver) const;

private:

    void setCircle(const Circle& circle, Scalar innerValue, ScalarFiniteVolumeField& field) const;
    void setCircle(const Circle& circle, const Vector2D& innerValue, VectorFiniteVolumeField& field) const;

    void setRotating(const std::string& function, Scalar amplitude, const Vector2D& center, ScalarFiniteVolumeField& field) const;
    void setRotating(const std::string& xFunction, const std::string& yFunction, const Vector2D& amplitude, const Vector2D& center, VectorFiniteVolumeField& field) const;

    boost::property_tree::ptree caseInput_;
    boost::property_tree::ptree boundaryInput_;
    boost::property_tree::ptree initialConditionInput_;

};

#endif
