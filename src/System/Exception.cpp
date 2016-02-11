#include "Exception.h"

Exception::Exception(const std::string &className, const std::string &methodName, const std::string &description)
    :
      message_(className + "::" + methodName + " -> " + description)
{

}
