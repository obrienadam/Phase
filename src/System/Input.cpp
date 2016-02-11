#include <boost/property_tree/info_parser.hpp>

#include "Input.h"

Input::Input(const std::string &caseFileName, const std::string &outputPath)
    :
      caseFileName(caseFileName),
      outputPath(outputPath)
{

}

void Input::parseInputFile()
{
    using namespace boost::property_tree;

    read_info(caseFileName, *((ptree*)this));
}
