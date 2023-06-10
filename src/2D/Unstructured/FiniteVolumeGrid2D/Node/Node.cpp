#include "FiniteVolumeGrid2D/FiniteVolumeGrid2D.h"

#include "Node.h"

Node::Node(Scalar x, Scalar y, const FiniteVolumeGrid2D &grid)
    : Point2D(x, y), grid_(grid), id_(grid.nodes().size()) {}

Node::Node(const Point2D &point, const FiniteVolumeGrid2D &grid)
    : Node(point.x, point.y, grid) {}

bool Node::isBoundaryNode() const {
  for (const Cell &c : cells_)
    for (const auto &nb : c.neighbours())
      if ((nb.face().lNode().id() == id_ || nb.face().rNode().id() == id_) &&
          nb.face().isBoundary())
        return true;

  return false;
}

std::vector<Scalar> Node::volumeWeights() const {
  Scalar sumW = 0.;
  std::vector<Scalar> weights;
  weights.reserve(4);

  for (const Cell &cell : cells()) {
    Scalar w = 1. / cell.volume();
    sumW += w;
    weights.push_back(w);
  }

  for (Scalar &w : weights)
    w /= sumW;

  return weights;
}

std::vector<Scalar> Node::distanceWeights() const {
  Scalar sumW = 0.;
  std::vector<Scalar> weights;

  for (const Cell &cell : cells()) {
    Scalar w = 1. / (*this - cell.centroid()).mag();
    sumW += w;
    weights.push_back(w);
  }

  for (Scalar &w : weights)
    w /= sumW;

  return weights;
}

bool Node::cellBounded(const Point2D &pt) const {
  std::vector<Point2D> centroids;
  for (const Cell &cell : cells())
    centroids.push_back(cell.centroid());

  return Polygon::convexHull(centroids.begin(), centroids.end()).isInside(pt);
}
