#include "CommandLine.h"

CommandLine::CommandLine()
    :
      desc_("Options")
{
    using namespace std;
    using namespace boost::program_options;

    desc_.add_options()
            ("help,h", "Print help message")
            ("case,c", value<string>()->value_name("<file>"), "Specify case file")
            ("output,o", value<string>()->value_name("<dir>"), "Specify output path");
}

CommandLine::CommandLine(int argc, const char *argv[], Input &input)
    :
      CommandLine()
{
    parseArguments(argc, argv, input);
}

void CommandLine::parseArguments(int argc, const char *argv[], Input& input)
{
    using namespace std;
    using namespace boost::program_options;

    store(parse_command_line(argc, argv, desc_), vm_);

    if(vm_.count("help"))
    {
        cout << "Phase: A multiphase 2D flow solver for problems involving immersed boundary methods.\n"
             << "Usage: phase [options]...\n"
             << desc_ << "\n";
    }

    if(vm_.count("case"))
    {
        input.caseDirectory = vm_["case"].as<string>();
    }

    if(vm_.count("output"))
    {
        input.outputPath = vm_["output"].as<string>();
    }
}
