#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>

#include "Types.h"

template <class T>
class Field : public std::vector<T>
{
public:

    Field(size_t size = 0, const T& initialValue = T(), const std::string& name = "N/A");
    Field(const Field<T>& other) : std::vector<T>(other), name_(other.name_) {}

    const std::string& name() const { return name_; }

protected:

    std::string name_;

};

#include "Field.tpp"

#endif
