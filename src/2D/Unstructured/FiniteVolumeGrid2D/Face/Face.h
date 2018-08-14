#ifndef PHASE_FACE_H
#define PHASE_FACE_H

#include "../Node/Node.h"
#include "../Cell/Cell.h"

class FiniteVolumeGrid2D;

class Face
{
public:

    enum Type
    {
        INTERIOR, BOUNDARY
    };

    Face(Label lNodeId, Label rNodeId, const FiniteVolumeGrid2D &grid, Type type = INTERIOR);

    void setType(Type type)
    { type_ = type; }

    Type type() const
    { return type_; }

    const Point2D &centroid() const
    { return centroid_; }

    const Vector2D &norm() const
    { return normal_; }

    Vector2D polarOutwardNorm(const Point2D& point) const;

    const Vector2D &tan() const
    { return tangent_; }

    Vector2D outwardNorm(const Point2D &point) const;

    Vector2D outwardNorm() const;

    Label id() const
    { return id_; }

    void setId(Label id)
    { id_ = id; }

    bool isInterior() const
    { return type_ == INTERIOR; }

    bool isBoundary() const
    { return type_ == BOUNDARY; }

    const Node &lNode() const;

    const Node &rNode() const;

    const Cell &lCell() const;

    const Cell &rCell() const;

    Scalar volumeWeight() const;

    Scalar distanceWeight() const;

    std::vector<Ref<const Cell>> cells() const;

    void addCell(const Cell &cell);

    std::string info() const;

    //- Grid
    const FiniteVolumeGrid2D &grid() const
    { return grid_; }

    //- Interpolation

private:
    Type type_;

    Point2D centroid_;

    Vector2D normal_, tangent_;

    Label id_;

    std::pair<Label, Label> nodeIds_;

    std::vector<Label> cellIds_;

    const FiniteVolumeGrid2D &grid_;
};

#endif
