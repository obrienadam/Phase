#include "Discretizer.h"

Term operator+(const Term &lhs, const Term &rhs)
{
    Term result;

    result.preallocate(lhs.triplets_.size() + rhs.triplets_.size());

    for(const Triplet& triplet: lhs.triplets_)
        result.triplets_.push_back(triplet);

    for(const Triplet& triplet: rhs.triplets_)
        result.triplets_.push_back(triplet);


    return result;
}

void Term::preallocate(Size size)
{
    triplets_.reserve(size);
    boundaries_ = SparseVector::Zero(size);
    sources_ = SparseVector::Zero(size);
}
