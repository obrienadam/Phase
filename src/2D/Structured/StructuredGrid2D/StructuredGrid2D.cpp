#include "StructuredGrid2D.h"

StructuredGrid2D::StructuredGrid2D(Size nNodesI, Size nNodesJ, Scalar lx, Scalar ly)
{
    init(nNodesI, nNodesJ, lx, ly);
}

void StructuredGrid2D::init(Size nNodesI, Size nNodesJ, Scalar lx, Scalar ly)
{
    nNodesI_ = nNodesI;
    nNodesJ_ = nNodesJ;
    lx_ = lx;
    ly_ = ly;

    Scalar dx = lx_ / (nNodesI_ - 1);
    Scalar dy = ly_ / (nNodesJ_ - 1);

    for (auto j = 0; j < nNodesJ_; ++j)
        for (auto i = 0; i < nNodesI_; ++i)
            nodes_.push_back(Node(Point2D(i * dx, j * dy), i, j, nodes_.size()));

    ownership_.assign(nodes_.size(), _comm.rank());
    globalIds_.resize(nodes_.size());
    std::iota(globalIds_.begin(), globalIds_.end(), 0);
}

Scalar StructuredGrid2D::dxe(size_t i, size_t j) const
{
    return node(i + 1, j).x - node(i, j).x;
}

Scalar StructuredGrid2D::dxw(size_t i, size_t j) const
{
    return node(i - 1, j).x - node(i, j).x;
}

Scalar StructuredGrid2D::dxn(size_t i, size_t j) const
{
    return node(i, j + 1).y - node(i, j).y;
}

Scalar StructuredGrid2D::dxs(size_t i, size_t j) const
{
    return node(i, j - 1).y - node(i, j).y;
}