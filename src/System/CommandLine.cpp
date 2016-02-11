#include "CommandLine.h"

CommandLine::CommandLine()
    :
      desc_("Options")
{
    using namespace boost::program_options;

    desc_.add_options()
            ("help,h", "Print help message")
            ("case,c", value<std::string>()->value_name("<file>"), "Specify case file")
            ("output,o", value<std::string>()->value_name("<dir>"), "Specify output path");
}

CommandLine::CommandLine(int argc, const char *argv[], Input &input)
    :
      CommandLine()
{
    parseArguments(argc, argv, input);
}

void CommandLine::parseArguments(int argc, const char *argv[], Input& input)
{
    using namespace boost::program_options;

    try
    {
        store(parse_command_line(argc, argv, desc_), vm_);

        if(vm_.count("help"))
        {
            std::cout << "Phase: A multiphase 2D flow solver for problems involving immersed boundary methods.\n"
                      << "Usage: phase [options]...\n"
                      << desc_ << "\n";
        }

        if(vm_.count("case"))
        {
            input.caseFileName = vm_["case"].as<std::string>();
        }

        if(vm_.count("output"))
        {
            input.outputPath = vm_["output"].as<std::string>();
        }
    }
    catch (const error& e)
    {
        std::cerr << "phase: fatal error: " << e.what() << "\n";
        exit(-1);
    }
}
