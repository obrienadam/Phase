#include "FiniteVolumeGrid2D.h"
#include "Exception.h"

FiniteVolumeGrid2D::FiniteVolumeGrid2D(size_t nNodes, size_t nCells, size_t nFaces)
{
    nodes_.reserve(nNodes);
    cells_.reserve(nCells);
    faces_.reserve(nFaces);
    patches_.reserve(4);
}

size_t FiniteVolumeGrid2D::createFace(size_t lNodeId, size_t rNodeId, Face::Type type)
{
    using namespace std;
    typedef pair<size_t, size_t> Key;

    Face face(lNodeId, rNodeId, nodes_, type);
    face.id_ = faces_.size();

    faces_.push_back(face);

    Key key1(lNodeId, rNodeId), key2(rNodeId, lNodeId);

    faceDirectory_.insert(pair<Key, size_t>(key1, face.id()));
    faceDirectory_.insert(pair<Key, size_t>(key2, face.id()));

    if(type == Face::INTERIOR)
        interiorFaces_.push_back(Ref<const Face>(faces_.back()));
    else if (type == Face::BOUNDARY)
        boundaryFaces_.push_back(Ref<const Face>(faces_.back()));

    return face.id();
}

size_t FiniteVolumeGrid2D::createCell(const std::vector<size_t> &faceIds)
{
    Cell cell(faceIds, faces_);
    cell.id_ = cells_.size();
    cell.globalIndex_ = cell.id();
    cells_.push_back(cell);

    for(size_t faceId: faceIds)
        faces_[faceId].addCell(cells_.back());

    return cell.id();
}

size_t FiniteVolumeGrid2D::findFace(size_t lNodeId, size_t rNodeId) const
{
    return faceDirectory_.find(std::pair<size_t, size_t>(lNodeId, rNodeId))->second;
}

std::string FiniteVolumeGrid2D::gridInfo() const
{
    using namespace std;

    return "Finite volume grid info:\n"
            "------------------------\n"
            "Number of nodes: " + to_string(nodes_.size()) + "\n"
            "Number of cells: " + to_string(cells_.size()) + "\n"
            "Number of interior faces: " + to_string(interiorFaces_.size()) + "\n"
            "Number of boundary faces: " + to_string(boundaryFaces_.size()) + "\n"
            "Number of faces total: " + to_string(faces_.size()) + "\n"
            "Bounding box: " + bBox_.toString() + "\n";
}

void FiniteVolumeGrid2D::computeCellAdjacency()
{
    for(Cell& cell: cells_)
        cell.computeCellAdjacency();
}

size_t FiniteVolumeGrid2D::computeGlobalIndices()
{
    size_t globalIndex = 0;

    for(Cell& cell: cells_)
        cell.globalIndex_ = cell.isActive() ? globalIndex++ : Cell::INACTIVE;

    return globalIndex;
}

const Cell& FiniteVolumeGrid2D::findContainingCell(const Point2D &point, const Cell &guess) const
{
    if(guess.isInCell(point))
        return guess;

    for(const InteriorLink &nb: guess.neighbours())
    {
        if(nb.cell().isInCell(point))
            return nb.cell();
    }

    for(const InteriorLink &nb: guess.neighbours())
        return findContainingCell(point, nb.cell());

    throw Exception("FiniteVolumeGrid2D", "findContainingCell", "specified point is not contained within any cell.");
}

//- Protected methods

void FiniteVolumeGrid2D::computeBoundingBox()
{
    bBox_ = BoundingBox(nodes_.data(), nodes_.size());
}

void FiniteVolumeGrid2D::applyPatch(const std::string &patchName, const std::vector<Ref<Face> > &faces)
{
    patches_.push_back(Patch(patches_.size(), patchName));

    for(Face &face: faces)
        patches_.back().addFace(face);
}
