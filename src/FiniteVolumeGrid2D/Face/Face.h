#ifndef FACE_H
#define FACE_H

#include <vector>
#include <memory>

#include "Node.h"
#include "Patch.h"

class Cell;

class Face
{
public:

   enum Type{INTERIOR, BOUNDARY};

   Face(size_t lNodeId, size_t rNodeId, const std::vector<Node>& nodes, Type type = INTERIOR);

   Type type() const { return type_; }

   const Point2D& centroid() const { return centroid_; }
   const Vector2D& norm() const { return normal_; }
   const Vector2D& tan() const { return tangent_; }
   Vector2D outwardNorm(const Point2D &point) const;

   const Patch &patch() const { return *patchPtr_; }

   Label id() const { return id_; }

   bool isInterior() const { return type_ == INTERIOR; }
   bool isBoundary() const { return type_ == BOUNDARY; }

   const Node& lNode() const { return nodes_[0]; }
   const Node& rNode() const { return nodes_[1]; }

   const Cell& lCell() const { return cells_[0]; }
   const Cell& rCell() const { return cells_[1]; }

   void addCell(const Cell& cell);

   std::string info() const;

private:

    void addToPatch(const Patch &patch) const;
    void changePatch(const Patch &patch) const;

    Type type_;

    Point2D centroid_;
    Vector2D normal_, tangent_;

    mutable const Patch *patchPtr_;

    Label id_;

    std::vector< Ref<const Cell> > cells_;
    std::vector< Ref<const Node> > nodes_;

    friend class FiniteVolumeGrid2D;
    friend class Patch;
};

#endif
