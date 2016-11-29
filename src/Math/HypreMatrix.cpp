#include <HYPRE.h>

#include "HypreMatrix.h"
#include "Exception.h"

HypreMatrix::HypreMatrix(const Communicator &comm,
                         int lOwnershipRange,
                         int uOwnershipRange)
{
    HYPRE_IJMatrixCreate(comm.comm(), lOwnershipRange, uOwnershipRange, lOwnershipRange, uOwnershipRange, &mat_);
    HYPRE_IJMatrixSetObjectType(mat_, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(mat_);
}

HypreMatrix::~HypreMatrix()
{
    HYPRE_IJMatrixDestroy(mat_);
}

void HypreMatrix::addValues(std::vector<int> &rows, std::vector<int> &nCols, std::vector<int> &cols, std::vector<double> &vals)
{
    if(!rows.size() == nCols.size())
        throw Exception("HypreMatrix", "assemble", "rows.size() != nCols.size().");
    else if(cols.size() != vals.size())
        throw Exception("HypreMatrix", "assemble", "cols.size() != vals.size().");

    HYPRE_IJMatrixAddToValues(mat_, rows.size(), nCols.data(), rows.data(), cols.data(), vals.data());
}

HypreMatrix &HypreMatrix::operator+=(const HypreMatrix &rhs)
{

}
