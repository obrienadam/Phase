#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <vector>
#include <map>

#include "Input.h"

class CommandLine
{
public:

    CommandLine();
    CommandLine(int argc, char* argv[]);

    void setOptions(const std::map<std::string, std::string>& options);
    void parseArguments(int argc, char *argv[]);

    std::string getOption(const std::string& option);

    int argc() const { return argc_; }
    char** argv() const { return argv_; }

private:

    void printHelpMessage(const char programName[]);
    void printVersrionInfo(const char programName[]);

    int argc_;
    char **argv_;

    std::map<std::string, std::string> options_;
    std::map<std::string, std::string> parsedArgs_;

};

#endif
