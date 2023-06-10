#ifndef PHASE_EXCEPTION_H
#define PHASE_EXCEPTION_H

#include <exception>
#include <string>

class Exception : public std::exception {
public:
  explicit Exception(const std::string &className,
                     const std::string &methodName,
                     const std::string &description);
  virtual ~Exception() throw() {}

  virtual const char *what() const throw() { return message_.c_str(); }

protected:
  std::string message_;
};

#endif
