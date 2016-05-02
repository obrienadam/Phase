#include "ScalarFiniteVolumeField.h"

template<>
ScalarFiniteVolumeField& ScalarFiniteVolumeField::operator =(const SparseVector& rhs)
{
    auto &self = *this;

    for(const Cell &cell: self.grid.activeCells())
        self[cell.id()] = rhs[cell.globalIndex()];

    return self;
}

//- Protected methods

template<>
void ScalarFiniteVolumeField::setBoundaryRefValues(const Input &input)
{
    using namespace std;

    auto &self = *this;

    //- Check for a default patch

    auto defBoundary = input.boundaryInput().get_child_optional("Boundaries." + Field<Scalar>::name + ".*");

    if(defBoundary)
    {
        for(int i = 0; i < grid.patches().size(); ++i)
           boundaryRefValues_.push_back(input.boundaryInput().get<Scalar>("Boundaries." + Field<Scalar>::name + ".*.value"));
    }
    else
    {
        //- Boundary condition association to patches is done by patch id
        for(const Patch& patch: grid.patches())
        {
            string root = "Boundaries." + Field<Scalar>::name + "." + patch.name;
            boundaryRefValues_.push_back(input.boundaryInput().get<Scalar>(root + ".value"));
        }
    }

    for(const Patch& patch: grid.patches())
    {
        for(const Face& face: patch.faces())
            self.faces()[face.id()] = boundaryRefValue(face.id());
    }
}

//- External functions

ScalarFiniteVolumeField operator*(const ScalarFiniteVolumeField& lhs, ScalarFiniteVolumeField rhs)
{
    rhs *= lhs;
    return rhs;
}
