#include <iostream>

#include "CommandLine.h"
#include "Exception.h"

CommandLine::CommandLine()
{
    using namespace std;

    options_ = map<string, string>{
        {"--help", "Displays this help message"},
        {"--version", "Displays version information"}
    };
}

CommandLine::CommandLine(int argc, char *argv[])
    :
      CommandLine()
{
    parseArguments(argc, argv);
}

void CommandLine::setOptions(const std::map<std::string, std::string> &options)
{
    options_ = options;
}

void CommandLine::parseArguments(int argc, char *argv[])
{
    argc_ = argc;
    argv_ = argv;

    for(int argNo = 1; argNo < argc; argNo += 2)
    {
        auto it = options_.find(argv[argNo]);

        if(it == options_.end())
            throw Exception("CommandLine", "parseArguments", "no such option \"" + std::string(argv[argNo]) + "\".");
        else if(it->first == "--help")
            printHelpMessage(argv[0]);
        else if(it->first == "--version")
            printVersrionInfo(argv[0]);
        else
            parsedArgs_[argv[argNo]] = argv[argNo + 1];
    }
}

std::string CommandLine::getOption(const std::string &option)
{
    auto it = parsedArgs_.find(option);

    if(it == parsedArgs_.end())
        throw Exception("CommandLine", "getOption", "unrecognized option \"" + option + "\".");

    return it->second;
}

//- Private

void CommandLine::printHelpMessage(const char programName[])
{
    printf("Usage: %s [options]\n", programName);
    printf("Must be executed from the top level of a case directory.\n");
    printf("Options:\n");

    for(const auto &opt: options_)
        printf("  %-23s%s\n", opt.first.c_str(), opt.second.c_str());

    exit(0);
}

void CommandLine::printVersrionInfo(const char programName[])
{
    printf("%s 0.0\n", programName);
    printf("Copyright (C) 2016 Adam R. O'Brien, a.obrien@mail.utoronto.ca\n");
    printf("This is free software; see the source for copying conditions.  There is NO\n");
    printf("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");

    exit(0);
}
