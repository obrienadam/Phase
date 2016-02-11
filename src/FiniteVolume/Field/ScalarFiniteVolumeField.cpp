#include "ScalarFiniteVolumeField.h"

ScalarFiniteVolumeField& ScalarFiniteVolumeField::operator =(const SparseVector& rhs)
{
    for(int j = 0, endJ = grid_.nCellsJ(); j < endJ; ++j)
        for(int i = 0, endI = grid_.nCellsI(); i < endI; ++i)
        {
            int row = grid_.cellIndex(i, j);

            if(row != FiniteVolumeGrid2D::INACTIVE)
                operator ()(i, j) = rhs[row];
        }

    return *this;
}
