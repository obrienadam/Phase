#ifndef PHASE_FIELD_H
#define PHASE_FIELD_H

#include <string>
#include <vector>

#include "Types/Types.h"

template <class T> class Field : public std::vector<T> {
public:
  Field(size_t size = 0, const T &initialValue = T(),
        const std::string &name = "N/A");

  Field(const Field<T> &other) : std::vector<T>(other), name_(other.name_) {}

  const std::string &name() const { return name_; }

protected:
  std::string name_;
};

#include "Field.tpp"

#endif
