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

    virtual Field<T>& operator+=(const Field<T>& rhs);
    virtual Field<T>& operator-=(const Field<T>& rhs);

    template <class T2>
    Field<T>& operator*=(const Field<T2>& rhs);

    std::string name;

protected:

    std::vector<T> prev_;

};

template <class T>
Field<T> operator+(Field<T> lhs, const Field<T>& rhs);

template <class T>
Field<T> operator-(Field<T> lhs, const Field<T>& rhs);

#include "Field.tpp"

#endif
