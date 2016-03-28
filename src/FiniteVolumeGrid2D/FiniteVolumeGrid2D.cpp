#include "FiniteVolumeGrid2D.h"

FiniteVolumeGrid2D::FiniteVolumeGrid2D(size_t nNodes, size_t nCells, size_t nFaces)
{
    nodes.reserve(nNodes);
    cells.reserve(nCells);
    faces.reserve(nFaces);
    patches_.reserve(4);
}

size_t FiniteVolumeGrid2D::createFace(size_t lNodeId, size_t rNodeId, Face::Type type)
{
    using namespace std;
    typedef pair<size_t, size_t> Key;

    Face face(lNodeId, rNodeId, nodes, type);
    face.id_ = faces.size();

    faces.push_back(face);

    Key key1(lNodeId, rNodeId), key2(rNodeId, lNodeId);

    faceDirectory_.insert(pair<Key, size_t>(key1, face.id()));
    faceDirectory_.insert(pair<Key, size_t>(key2, face.id()));

    if(type == Face::INTERIOR)
        interiorFaces_.push_back(face.id());
    else if (type == Face::BOUNDARY)
        boundaryFaces_.push_back(face.id());

    return face.id();
}

size_t FiniteVolumeGrid2D::createCell(const std::vector<size_t> &faceIds)
{
    Cell cell(faceIds, faces);
    cell.id_ = cells.size();
    cell.globalIndex_ = cell.id();
    cells.push_back(cell);

    for(size_t faceId: faceIds)
        faces[faceId].addCell(cells.back());

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
            "Number of nodes: " + to_string(nodes.size()) + "\n"
            "Number of cells: " + to_string(cells.size()) + "\n"
            "Number of interior faces: " + to_string(interiorFaces_.size()) + "\n"
            "Number of boundary faces: " + to_string(boundaryFaces_.size()) + "\n"
            "Number of faces total: " + to_string(faces.size()) + "\n"
            "Bounding box: " + bBox_.toString() + "\n";
}

void FiniteVolumeGrid2D::computeCellAdjacency()
{
    for(Cell& cell: cells)
        cell.computeCellAdjacency();
}

ScalarFiniteVolumeField& FiniteVolumeGrid2D::addScalarField(const std::string &fieldName) const
{
    typedef std::pair< std::string, ScalarFiniteVolumeField> Key;

    return (scalarFields_.insert(Key(fieldName, ScalarFiniteVolumeField(*this, fieldName))).first)->second;
}

ScalarFiniteVolumeField& FiniteVolumeGrid2D::addScalarField(const Input& input, const std::string &fieldName) const
{
    typedef std::pair< std::string, ScalarFiniteVolumeField> Key;

    return (scalarFields_.insert(Key(fieldName, ScalarFiniteVolumeField(input, *this, fieldName))).first)->second;
}

VectorFiniteVolumeField& FiniteVolumeGrid2D::addVectorField(const Input &input, const std::string &fieldName) const
{
    typedef std::pair< std::string, VectorFiniteVolumeField> Key;

    return (vectorFields_.insert(Key(fieldName, VectorFiniteVolumeField(input, *this, fieldName))).first)->second;
}

VectorFiniteVolumeField& FiniteVolumeGrid2D::addVectorField(const std::string &fieldName) const
{
    typedef std::pair< std::string, VectorFiniteVolumeField> Key;

    return (vectorFields_.insert(Key(fieldName, VectorFiniteVolumeField(*this, fieldName))).first)->second;
}

//- Protected methods

void FiniteVolumeGrid2D::computeBoundingBox()
{
    bBox_ = BoundingBox(nodes.data(), nodes.size());
}

void FiniteVolumeGrid2D::applyPatch(const std::string &patchName, const std::vector<Ref<Face> > &faces)
{
    patches_.push_back(Patch(patches_.size(), patchName));

    for(Face &face: faces)
        patches_.back().addFace(face);
}
