#ifndef HYPRE_MATRIX
#define HYPRE_MATRIX

#include <vector>

#include <HYPRE_IJ_mv.h>
#include "Communicator.h"

class HypreMatrix
{
public:
    HypreMatrix(const Communicator& comm,
                int lOwnershipRange,
                int uOwnershipRange);

    ~HypreMatrix();

    void addValues(std::vector<int>& rows, std::vector<int>& nCols, std::vector<int>& cols, std::vector<double>& vals);

    HypreMatrix& operator+=(const HypreMatrix& rhs);

private:
    HYPRE_IJMatrix mat_;
    int numRows_, numCols_;
};

#endif
