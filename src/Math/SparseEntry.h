#ifndef PHASE_SPARSE_ENTRY_H
#define PHASE_SPARSE_ENTRY_H

#include "Types/Types.h"

class SparseEntry
{
public:

    SparseEntry() {}

    SparseEntry(Index row, Index col, Scalar val)
        : row(row), col(col), val(val)
    {}

    Index row, col;

    Scalar val;
};

#endif
