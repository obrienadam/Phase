#include "StructuredGrid2D.h"

StructuredGrid2D::StructuredGrid2D(size_t nNodesI, size_t nNodesJ, Scalar lx, Scalar ly)
        :
        nNodesI_(nNodesI),
        nNodesJ_(nNodesJ),
        lx_(lx),
        ly_(ly)
{
    Scalar dx = lx_ / (nNodesI_ - 1);
    Scalar dy = ly_ / (nNodesJ_ - 1);

    for(auto j = 0; j < nNodesJ_; ++j)
        for(auto i = 0; i < nNodesI_; ++i)
            nodes_.push_back(Point2D(i * dx, j * dy));
}