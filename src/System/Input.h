#ifndef INPUT_H
#define INPUT_H

#include <string>

#include <boost/property_tree/ptree.hpp>

class Input : public boost::property_tree::ptree
{
public:

    Input(const std::string& caseFileName = "case/case.info", const std::string& outputPath = "solution");

    void parseInputFile();

    std::string caseFileName, outputPath;

private:

};

#endif
