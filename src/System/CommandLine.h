#ifndef PHASE_COMMAND_LINE_H
#define PHASE_COMMAND_LINE_H

#include <unordered_map>

#include <boost/program_options.hpp>

class CommandLine
{
public:

    CommandLine();

    CommandLine(int argc, char* argv[]);

    void addSwitch(const std::string& opt, const std::string& desc);

    template<class T>
    void add(const std::string& opt, const std::string& desc, bool optional = true)
    {
        if(optional)
            desc_.add_options()(opt.c_str(), boost::program_options::value<T>(), desc.c_str());
        else
            desc_.add_options()(opt.c_str(), boost::program_options::value<T>()->required(), desc.c_str());
    }

    template<class T>
    const T& get(const std::string& opt) const
    {
        return vm_[opt].as<T>();
    }

    boost::program_options::options_description_easy_init addOptions();

    void parseArguments(int argc, char *argv[]);

    std::string getOption(const std::string& option);


private:

    boost::program_options::options_description desc_;

    boost::program_options::variables_map vm_;

};

#endif
