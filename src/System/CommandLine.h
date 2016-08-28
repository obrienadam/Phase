#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <vector>
#include <map>

#include "Input.h"

class CommandLine
{
public:

    CommandLine();
    CommandLine(int argc, const char* argv[], Input& input);

    void parseArguments(int argc, const char* argv[], Input &input);

    int argc() const { return argc_; }
    const char** argv() const { return argv_; }

private:

    void printHelpMessage(const char programName[]);
    void printVersrionInfo(const char programName[]);

    int argc_;
    const char **argv_;

    std::map<std::string, std::string> options_;
    std::map<std::string, std::string> parsedArgs_;

};

#endif
