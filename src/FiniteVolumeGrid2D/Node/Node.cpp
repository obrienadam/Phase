#include "Node.h"
#include "FiniteVolumeGrid2D.h"

Node::Node(Scalar x, Scalar y, const FiniteVolumeGrid2D &grid)
        :
        Point2D(x, y),
        cells_(grid.cells())
{
    id_ = grid.nodes().size();
}

Node::Node(const Point2D &point, const FiniteVolumeGrid2D &grid)
        :
        Node(point.x, point.y, grid)
{

}

void Node::addCell(const Cell &cell)
{
    cellIds_.push_back(cell.id());
}

const std::vector<Ref<const Cell> > Node::cells() const
{
    using namespace std;
    vector<Ref<const Cell>> cells;

    for (Label id: cellIds_)
        cells.push_back(cref(cells_[id]));

    return cells;
}
