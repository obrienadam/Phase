#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>
#include <memory>

#include "Types.h"

template <class T>
class Field : public std::vector<T>
{
public:

    Field(size_t size = 0, const T& initialValue = T(), const std::string& name = "N/A");
    Field(const Field<T>& other) : std::vector<T>(other), name(other.name), prevFieldPtr_(nullptr) {}

    void save();

    Field<T>& prev() { return *prevFieldPtr_; }
    Field<T>& prev() const { return *prevFieldPtr_; }
    bool previousFieldAvailable() const { return prevFieldPtr_ != nullptr; }

    virtual Field<T>& operator+=(const Field<T>& rhs);
    virtual Field<T>& operator-=(const Field<T>& rhs);

    template <class T2>
    Field<T>& operator*=(const Field<T2>& rhs);

    T max() const;

    std::string name;

protected:

    std::shared_ptr< Field<T> > prevFieldPtr_;

};

template <class T>
Field<T> operator+(Field<T> lhs, const Field<T>& rhs);

template <class T>
Field<T> operator-(Field<T> lhs, const Field<T>& rhs);

#include "Field.tpp"

#endif
