#ifndef PHASE_NOT_IMPLEMENTED_EXCEPTION_H
#define PHASE_NOT_IMPLEMENTED_EXCEPTION_H

#include "Exception.h"

class NotImplementedException : public Exception
{
public:
    explicit NotImplementedException(const std::string &className, const std::string &methodName)
            :
            Exception(className, methodName, "is not implemented.")
    {}
};


#endif
