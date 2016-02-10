#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>

template <class T>
class Field : private std::vector<T>
{
public:

    Field(int sizeI = 1, int sizeJ = 1, const T& initialValue = 0, const std::string& name = "N/A");
    Field(int sizeI, int sizeJ, const std::string &name);

    T& operator()(int i, int j);
    const T& operator()(int i, int j) const;

    void fill(const T& value);
    void resize(int sizeI, int sizeJ);

    int size() const { return std::vector<T>::size(); }
    inline int sizeI() const { return sizeI_; }
    inline int sizeJ() const { return sizeJ_; }

    std::string name;

private:

    int sizeI_, sizeJ_;

};

#include "Field.tpp"

#endif
