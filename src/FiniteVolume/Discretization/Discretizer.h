#ifndef DISCRETIZER_H
#define DISCRETIZER_H

#include <vector>
#include "Types.h"
#include "SparseMatrix.h"

class Term
{
public:

private:
    std::vector<Triplet> triplets_;
    SparseVector boundaries_;
    SparseVector sources_;

    friend Term operator+(const Term &lhs, const Term &rhs);
};

Term operator+(const Term &lhs, const Term &rhs);

#endif
