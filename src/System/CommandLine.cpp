#include <iostream>

#include <boost/program_options.hpp>

#include "CommandLine.h"
#include "Exception.h"

CommandLine::CommandLine() {
  desc_.add_options()("help,h", "output this help message");
}

CommandLine::CommandLine(int argc, char *argv[]) : CommandLine() {
  parseArguments(argc, argv);
}

void CommandLine::addSwitch(const std::string &opt, const std::string &desc) {
  using namespace boost::program_options;
  desc_.add_options()(opt.c_str(), bool_switch()->default_value(false),
                      desc.c_str());
}

boost::program_options::options_description_easy_init
CommandLine::addOptions() {
  return desc_.add_options();
}

void CommandLine::parseArguments(int argc, char *argv[]) {
  boost::program_options::store(
      boost::program_options::parse_command_line(argc, argv, desc_), vm_);
  boost::program_options::notify(vm_);

  if (vm_.count("help")) {
    std::cout << desc_ << "\n";
    exit(0);
  }
}
