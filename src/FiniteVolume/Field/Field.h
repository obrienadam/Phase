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

    void save();

    std::vector<T>& prev() { return prev_; }
    const std::vector<T>& prev() const { return prev_; }
    bool previousFieldAvailable() const;

    std::string name;

protected:

    std::vector<T> prev_;

};

#include "Field.tpp"

#endif
