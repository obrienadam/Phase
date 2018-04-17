#include "StructuredGrid2D.h"

StructuredGrid2D::StructuredGrid2D(size_t nNodesI, size_t nNodesJ, Scalar lx, Scalar ly)
        :
        nNodesI_(nNodesI),
        nNodesJ_(nNodesJ),
        nCellsI_(nNodesI_ - 1),
        nCellsJ_(nNodesJ_ - 1),
        lx_(lx),
        ly_(ly)
{
    Scalar dx = lx_ / nCellsI_;
    Scalar dy = ly_ / nCellsJ_;

    for (auto j = 0; j < nNodesJ_; ++j)
        for (auto i = 0; i < nNodesI_; ++i)
            nodes_.push_back(Point2D(i * dx, j * dy));
}

Scalar StructuredGrid2D::vol(size_t i, size_t j) const
{

}

Scalar StructuredGrid2D::dxe(size_t i, size_t j) const
{

}

Scalar StructuredGrid2D::dxw(size_t i, size_t j) const
{

}

Scalar StructuredGrid2D::dxn(size_t i, size_t j) const
{

}

Scalar StructuredGrid2D::dxs(size_t i, size_t j) const
{

}