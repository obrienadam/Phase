#include <iostream>

#include "CommandLine.h"

CommandLine::CommandLine()
    :
      desc_("Options")
{
    using namespace boost::program_options;

    desc_.add_options()
            ("help", "Print help message")
            ("case", value<std::string>(), "Specify case file")
            ("output", value<std::string>(), "Specify output path");
}

CommandLine::CommandLine(int argc, const char *argv[], Input &input)
    :
      CommandLine()
{
    //parseArguments(argc, argv, input);
}

void CommandLine::parseArguments(int argc, const char *argv[], Input& input)
{
    using namespace boost::program_options;

    store(parse_command_line(argc, argv, desc_), vm_);

    notify(vm_);

    if(vm_.count("help"))
    {
        std::cout << "Phase: A multiphase 2D flow solver for problems involving immersed boundary methods.\n"
                  << "Usage: phase [options]...\n"
                  << desc_ << "\n";
        exit(0);
    }
    else
    {
        if(vm_.count("case"))
            input.caseDirectory = vm_["case"].as<std::string>();
        if(vm_.count("output"))
            input.outputPath = vm_["output"].as<std::string>();
    }
}
