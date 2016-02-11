#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <boost/program_options.hpp>

#include "Input.h"

class CommandLine
{
public:

    CommandLine();
    CommandLine(int argc, const char* argv[], Input& input);

    void parseArguments(int argc, const char* argv[], Input &input);

private:

    boost::program_options::options_description desc_;
    boost::program_options::variables_map vm_;

};

#endif
