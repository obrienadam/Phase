#ifndef INPUT_H
#define INPUT_H

#include <string>

#include <boost/property_tree/ptree.hpp>

class Input
{
public:

    Input(const std::string& caseDirectory = "case", const std::string& outputPath = "solution");

    void parseInputFile();

    std::string caseDirectory, outputPath;

    const boost::property_tree::ptree& caseInput() const { return caseInput_; }
    const boost::property_tree::ptree& boundaryInput() const { return boundaryInput_; }
    const boost::property_tree::ptree& initialConditionInput() const { return initialConditionInput_; }

private:

    boost::property_tree::ptree caseInput_;
    boost::property_tree::ptree boundaryInput_;
    boost::property_tree::ptree initialConditionInput_;

};

#endif
