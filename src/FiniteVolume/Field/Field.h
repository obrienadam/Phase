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

    std::string name;

protected:

};

#include "Field.tpp"

#endif
