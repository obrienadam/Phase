#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <exception>

class Exception : public std::exception
{
public:

    explicit Exception(const std::string& className, const std::string& methodName, const std::string& description);
    virtual ~Exception() throw() {}

    virtual const char* what() const throw() { return message_.c_str(); }

protected:

    std::string message_;
};

#endif
