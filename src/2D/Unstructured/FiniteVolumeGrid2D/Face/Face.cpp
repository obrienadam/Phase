#include <algorithm>

#include "FiniteVolumeGrid2D/FiniteVolumeGrid2D.h"
#include "System/Exception.h"

#include "Face.h"
#include "FiniteVolume/Discretization/Axisymmetric.h"

Face::Face(Label lNodeId, Label rNodeId, const FiniteVolumeGrid2D &grid,
           Type type)
    : type_(type), grid_(grid),
      nodes_(grid.nodes()[lNodeId], grid.nodes()[rNodeId]) {
  centroid_ = 0.5 * (lNode() + rNode());
  tangent_ = rNode() - lNode();
  normal_ = tangent_.normalVec();
  id_ = grid.faces().size();
  cells_.reserve(2);
}

Vector2D Face::polarOutwardNorm(const Point2D &point) const {
  Scalar r1 = lNode().x;
  Scalar r2 = rNode().x;
  Scalar z1 = lNode().y;
  Scalar z2 = rNode().y;

  Vector2D sf((r1 + r2) / 2. * (z2 - z1), (r2 * r2 - r1 * r1) / 2.);

  return dot(sf, centroid_ - point) > 0. ? sf : -sf;
}

Vector2D Face::polarOutwardNorm(const Point2D &point,
                                const Vector2D &zaxis) const {
  Vector2D ds1 = axi::map(lNode(), zaxis);
  Vector2D ds2 = axi::map(rNode(), zaxis);

  Scalar r1 = ds1.x;
  Scalar r2 = ds2.x;
  Scalar z1 = ds1.y;
  Scalar z2 = ds2.y;

  Vector2D sf((r1 + r2) / 2. * (z2 - z1), (r2 * r2 - r1 * r1) / 2.);

  return dot(sf, axi::map(centroid_, zaxis) - axi::map(point, zaxis)) > 0.
             ? sf
             : -sf;
}

Vector2D Face::outwardNorm(const Point2D &point) const {
  return dot(centroid_ - point, normal_) > 0. ? normal_ : -normal_;
}

Vector2D Face::outwardNorm() const { return outwardNorm(lCell().centroid()); }

Scalar Face::volumeWeight() const {
  Scalar v1 = rCell().volume();
  Scalar v2 = lCell().volume();
  return v1 / (v1 + v2);
}

Scalar Face::distanceWeight() const {
  Scalar l1 = (centroid_ - rCell().centroid()).mag();
  Scalar l2 = (centroid_ - lCell().centroid()).mag();
  return l1 / (l1 + l2);
}

void Face::addCell(const Cell &cell) {
  if (type_ == INTERIOR) {
    if (cells_.size() == 2)
      throw Exception(
          "Face", "addCell",
          "an interior face cannot be shared between more than two cells. " +
              info());
  } else if (type_ == BOUNDARY) {
    if (cells_.size() == 1)
      throw Exception(
          "Face", "addCell",
          "a boundary face cannot be shared by more than one cell. " + info());
  }

  cells_.emplace_back(cell);
}

std::string Face::info() const {
  using namespace std;

  return "Face info: "
         "Id: " +
         to_string(id_) + " " + "Node 1: " + to_string(lNode().id()) + ", " +
         "Node 2: " + to_string(rNode().id()) + ", " +
         "Centroid: " + centroid_.toString() + ", " +
         "Normal: " + normal_.toString();
}

//- Private methods
