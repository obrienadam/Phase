#include <fstream>

#include <boost/algorithm/string.hpp>

#include "System/Exception.h"

#include "FiniteVolume/Field/FiniteVolumeField.h"

//- Static

//- Parallel
template<class T>
std::vector<std::vector<T>> FiniteVolumeField<T>::sendBuffers_;

template<class T>
std::vector<std::vector<T>> FiniteVolumeField<T>::recvBuffers_;

//- Constructors

template<class T>
FiniteVolumeField<T>::FiniteVolumeField(const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                        const std::string &name,
                                        const T &val,
                                        bool faces,
                                        bool nodes,
                                        const std::shared_ptr<const CellGroup> &cellGroup,
                                        const std::shared_ptr<IndexMap> &indexMap)
    :
      Field<T>::Field(grid->cells().size(), val, name),
      grid_(grid),
      indexMap_(indexMap)
{
    cellGroup_ = cellGroup;

    if (faces)
        faces_.resize(grid_->faces().size(), val);

    if (nodes)
        nodes_.resize(grid_->nodes().size(), val);
}

template<class T>
FiniteVolumeField<T>::FiniteVolumeField(const Input &input,
                                        const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                        const std::string &name,
                                        const T &val,
                                        bool faces,
                                        bool nodes,
                                        const std::shared_ptr<const CellGroup> &cellGroup,
                                        const std::shared_ptr<IndexMap> &indexMap)
    :
      FiniteVolumeField(grid, name, val, faces, nodes, cellGroup, indexMap)
{
    setBoundaryTypes(input);
    setBoundaryRefValues(input);
}

//- Public methods

template<class T>
void FiniteVolumeField<T>::fill(const T &val)
{
    std::fill(this->begin(), this->end(), val);
    std::fill(faces_.begin(), faces_.end(), val);
    std::fill(nodes_.begin(), nodes_.end(), val);
}

template<class T>
void FiniteVolumeField<T>::fill(const T &val, const CellGroup &group)
{
    for (const Cell &cell: group)
        (*this)[cell.id()] = val;
}

template<class T>
void FiniteVolumeField<T>::fillInterior(const T &val)
{
    std::fill(this->begin(), this->end(), val);

    if (!faces_.empty())
    {
        for (const Face &face: grid_->interiorFaces())
            faces_[face.id()] = val;
    }
}

template<class T>
void FiniteVolumeField<T>::assign(const FiniteVolumeField<T> &field)
{
    this->assign(field.begin(), field.end());
    faces_.assign(field.faces_.begin(), field.faces_.end());
    nodes_.assign(field.nodes_.begin(), field.nodes_.end());
    patchBoundaries_ = field.patchBoundaries_;
    grid_ = field.grid_;
    cellGroup_ = field.cellGroup_;
}

template<class T>
void FiniteVolumeField<T>::copyBoundaryTypes(const FiniteVolumeField &other)
{
    patchBoundaries_ = other.patchBoundaries_;
}

template<class T>
typename FiniteVolumeField<T>::BoundaryType FiniteVolumeField<T>::boundaryType(const FaceGroup &patch) const
{
    auto it = patchBoundaries_.find(patch.name());
    return it == patchBoundaries_.end() ? NORMAL_GRADIENT : it->second.first;
}

template<class T>
typename FiniteVolumeField<T>::BoundaryType FiniteVolumeField<T>::boundaryType(const Face &face) const
{
    return boundaryType(grid_->patch(face));
}

template<class T>
T FiniteVolumeField<T>::boundaryRefValue(const FaceGroup &patch) const
{
    return patchBoundaries_.find(patch.name())->second.second;
}

template<class T>
T FiniteVolumeField<T>::boundaryRefValue(const Face &face) const
{
    return boundaryRefValue(grid_->patch(face));
}

template<class T>
template<class TFunc>
void FiniteVolumeField<T>::interpolateFaces(const TFunc &alpha)
{
    auto &self = *this;

    for (const Face &face: grid_->interiorFaces())
    {
        Scalar g = alpha(face);
        self(face) = g * self(face.lCell()) + (1. - g) * self(face.rCell());
    }

    setBoundaryFaces();
}

template<class T>
void FiniteVolumeField<T>::interpolateFaces(InterpolationType type)
{
    switch (type)
    {
    case VOLUME:
        interpolateFaces([](const Face &face)
        {
            return face.volumeWeight();
        });
        break;
    case DISTANCE:
        interpolateFaces([](const Face &face)
        {
            return face.distanceWeight();
        });
        break;
    }
}

template<class T>
void FiniteVolumeField<T>::setBoundaryFaces()
{
    auto &self = *this;

    for (const FaceGroup &patch: grid_->patches())
    {
        switch (boundaryType(patch))
        {
        case FIXED:
            break;
        case NORMAL_GRADIENT:
        case SYMMETRY:
            for (const Face &face: patch)
                faces_[face.id()] = self[face.lCell().id()];
            break;
        }
    }
}

template<class T>
void FiniteVolumeField<T>::setBoundaryFaces(BoundaryType bType, const std::function<T(const Face &face)> &fcn)
{
    auto &self = *this;
    for (const FaceGroup &patch: grid_->patches())
    {
        if (boundaryType(patch) == bType)
            for (const Face &face: patch)
                self(face) = fcn(face);
    }
}

template<class T>
void FiniteVolumeField<T>::interpolateNodes()
{
    for(const Node &node: grid_->nodes())
    {
        for(const Cell& cell: node.cells())
        {

        }
    }
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::savePreviousTimeStep(Scalar timeStep, int nPreviousFields)
{
    std::shared_ptr<FiniteVolumeField<T>> tmp;

    if(previousTimeSteps_.size() >= nPreviousFields)
    {
        tmp = previousTimeSteps_.back().second;
        *tmp = *this;
    }
    else
        tmp = std::make_shared<FiniteVolumeField<T>>(*this);

    tmp->clearHistory();

    previousTimeSteps_.emplace_front(timeStep, tmp);
    previousTimeSteps_.resize(nPreviousFields, previousTimeSteps_.back());

    return *previousTimeSteps_.front().second;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::savePreviousIteration()
{
    if(previousIteration_)
        *previousIteration_ = *this;
    else
        previousIteration_ = std::make_shared<FiniteVolumeField<T>>(*this);

    previousIteration_->clearHistory();
    return *previousIteration_;
}

template<class T>
void FiniteVolumeField<T>::clearHistory()
{
    previousIteration_ = nullptr;
    previousTimeSteps_.clear();
}

//- Parallel

template<class T>
void FiniteVolumeField<T>::sendMessages()
{
    sendBuffers_.resize(grid_->comm().nProcs());
    recvBuffers_.resize(grid_->comm().nProcs());

    //- Post recvs
    for(int proc = 0; proc < grid_->comm().nProcs(); ++proc)
    {
        if(proc == grid_->comm().rank())
            continue;

        recvBuffers_[proc].resize(grid_->bufferGroups()[proc].size());
        grid_->comm().irecv(proc, recvBuffers_[proc], proc);
    }

    //- Post sends
    for(int proc = 0; proc < grid_->comm().nProcs(); ++proc)
    {
        if(proc == grid_->comm().rank())
            continue;

        const CellSet &sendCells = grid_->sendGroups()[proc];

        sendBuffers_[proc].resize(sendCells.size());
        std::transform(sendCells.begin(), sendCells.end(), sendBuffers_[proc].begin(),
                       [this](const Cell &c) { return (*this)(c); });

        grid_->comm().isend(proc, sendBuffers_[proc], grid_->comm().rank());
    }

    //- Sync
    grid_->comm().waitAll();

    //- Unload recv buffers
    for(int proc = 0; proc < grid_->comm().nProcs(); ++proc)
    {
        if(proc == grid_->comm().rank())
            continue;

        const std::vector<T> &buff = recvBuffers_[proc];

        int i = 0;
        for(const Cell &c: grid_->bufferGroups()[proc])
            (*this)(c) = buff[i++];
    }
}

//- Operators

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator+=(const FiniteVolumeField &rhs)
{
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<T>());
    std::transform(faces_.begin(), faces_.end(), rhs.faces().begin(), faces_.begin(), std::plus<T>());
    std::transform(nodes_.begin(), nodes_.end(), rhs.nodes().begin(), nodes_.begin(), std::plus<T>());
    return *this;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator-=(const FiniteVolumeField &rhs)
{
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::minus<T>());
    std::transform(faces_.begin(), faces_.end(), rhs.faces().begin(), faces_.begin(), std::minus<T>());
    std::transform(nodes_.begin(), nodes_.end(), rhs.nodes().begin(), nodes_.begin(), std::minus<T>());
    return *this;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator*=(const FiniteVolumeField<Scalar> &rhs)
{
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), [](const T &lhs, Scalar rhs)
    {
        return lhs * rhs;
    });
    std::transform(faces_.begin(), faces_.end(), rhs.faces().begin(), faces_.begin(), [](const T &lhs, Scalar rhs)
    {
        return lhs * rhs;
    });
    std::transform(nodes_.begin(), nodes_.end(), rhs.nodes().begin(), nodes_.begin(), [](const T &lhs, const Scalar rhs)
    {
        return lhs * rhs;
    });

    return *this;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator/=(const FiniteVolumeField<Scalar> &rhs)
{
    std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), [](const T &lhs, Scalar rhs)
    {
        return lhs / rhs;
    });
    std::transform(faces_.begin(), faces_.end(), rhs.faces().begin(), faces_.begin(), [](const T &lhs, Scalar rhs)
    {
        return lhs / rhs;
    });
    std::transform(nodes_.begin(), nodes_.end(), rhs.nodes().begin(), nodes_.begin(), [](const T &lhs, Scalar rhs)
    {
        return lhs / rhs;
    });

    return *this;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator*=(Scalar rhs)
{
    std::transform(this->begin(), this->end(), this->begin(), [rhs](const T &lhs)
    {
        return lhs * rhs;
    });
    std::transform(faces_.begin(), faces_.end(), faces_.begin(), [rhs](const T &lhs)
    {
        return lhs * rhs;
    });
    std::transform(nodes_.begin(), nodes_.end(), nodes_.begin(), [rhs](const T &lhs)
    {
        return lhs * rhs;
    });

    return *this;
}

template<class T>
FiniteVolumeField<T> &FiniteVolumeField<T>::operator/=(Scalar rhs)
{
    std::transform(this->begin(), this->end(), this->begin(), [rhs](const T &lhs)
    {
        return lhs / rhs;
    });
    std::transform(faces_.begin(), faces_.end(), faces_.begin(), [rhs](const T &lhs)
    {
        return lhs / rhs;
    });
    std::transform(nodes_.begin(), nodes_.end(), nodes_.begin(), [rhs](const T &lhs)
    {
        return lhs / rhs;
    });

    return *this;
}

template<class T>
void FiniteVolumeField<T>::setGrid(const std::shared_ptr<const FiniteVolumeGrid2D> &grid)
{
    grid_ = grid;

    Field<T>::resize(grid_->cells().size());

    if (!faces_.empty())
        faces_.resize(grid_->faces().size());

    if (!nodes_.empty())
        nodes_.resize(grid_->nodes().size());

    cellGroup_ = nullptr;
}

//- Debug

template<class T>
void FiniteVolumeField<T>::writeToFile(const std::string &filename) const
{
    std::ofstream fout(filename);

    for (const T &val: *this)
        fout << val << "\n";

    fout.close();
}

//- Protected methods

template<class T>
void FiniteVolumeField<T>::setBoundaryTypes(const Input &input)
{
    using namespace std;


    std::string typeStr = input.boundaryInput().get<string>("Boundaries." + Field<T>::name() + ".*.type", "");
    BoundaryType boundaryType;
    std::vector<Scalar> coeffs;

    if (!typeStr.empty())
    {
        if (typeStr == "fixed")
            boundaryType = FIXED;
        else if (typeStr == "normal_gradient")
            boundaryType = NORMAL_GRADIENT;
        else if (typeStr == "symmetry")
            boundaryType = SYMMETRY;
        else if (typeStr == "outflow")
            boundaryType = OUTFLOW;
        else if (typeStr == "partial_slip")
            boundaryType = PARTIAL_SLIP;
        else
            throw Exception("FiniteVolumeField<T>", "setBoundaryTypes", "invalid boundary type \"" + typeStr + "\".");

        for (const FaceGroup &patch: grid_->patches())
        {
            patchBoundaries_[patch.name()] = std::make_pair(boundaryType, T());
        }
    }

    for (const FaceGroup &patch: grid_->patches())
    {
        typeStr = input.boundaryInput().get<std::string>(
                    "Boundaries." + Field<T>::name() + "." + patch.name() + ".type",
                    "");

        if (typeStr.empty())
            continue;

        if (typeStr == "fixed")
            boundaryType = FIXED;
        else if (typeStr == "normal_gradient")
            boundaryType = NORMAL_GRADIENT;
        else if (typeStr == "symmetry")
            boundaryType = SYMMETRY;
        else if (typeStr == "outflow")
            boundaryType = OUTFLOW;
        else
            throw Exception("FiniteVolumeField<T>", "setBoundaryTypes", "invalid boundary type \"" + typeStr + "\".");

        patchBoundaries_[patch.name()] = std::make_pair(boundaryType, T());
    }
}

//- External operators

template<class T>
FiniteVolumeField<T> operator+(FiniteVolumeField<T> lhs, const FiniteVolumeField<T> &rhs)
{
    lhs += rhs;
    return lhs;
}

template<class T>
FiniteVolumeField<T> operator-(FiniteVolumeField<T> lhs, const FiniteVolumeField<T> &rhs)
{
    lhs -= rhs;
    return lhs;
}

template<class T>
FiniteVolumeField<T> operator*(FiniteVolumeField<T> lhs, Scalar rhs)
{
    lhs *= rhs;
    return lhs;
}

template<class T>
FiniteVolumeField<T> operator*(Scalar lhs, FiniteVolumeField<T> rhs)
{
    rhs *= lhs;
    return rhs;
}

template<class T>
FiniteVolumeField<T> operator/(FiniteVolumeField<T> lhs, Scalar rhs)
{
    lhs /= rhs;
    return lhs;
}
