#include <functional>

#include <boost/property_tree/info_parser.hpp>

#include "Input.h"
#include "Solver.h"
#include "Exception.h"

Input::Input(const std::string &caseDirectory, const std::string &outputPath)
    :
      caseDirectory(caseDirectory),
      outputPath(outputPath)
{

}

void Input::parseInputFile()
{
    using namespace boost::property_tree;

    read_info(caseDirectory + "/case.info", caseInput_);
    read_info(caseDirectory + "/boundaries.info", boundaryInput_);
    read_info(caseDirectory + "/initialConditions.info", initialConditionInput_);
}

//- Private methods
