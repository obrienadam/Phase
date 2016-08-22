#include <fstream>

#include "ImmersedBoundary.h"
#include "Solver.h"
#include "GhostCellImmersedBoundary.h"

ImmersedBoundary::ImmersedBoundary(const Input &input, Solver &solver)
    :
      solver_(solver),
      cellStatus_(solver.addScalarField("cell_status"))
{
    try //- Lazy way to check if any immersed boundary input is present
    {
        input.boundaryInput().get_child("ImmersedBoundaries");
    }
    catch(...)
    {
        printf("No immersed boundaries present.\n");
        return;
    }

    for(const auto& ibObject: input.boundaryInput().get_child("ImmersedBoundaries"))
    {
        printf("Initializing immersed boundary object \"%s\".\n", ibObject.first.c_str());

        //- Initialize the geometry
        const std::string type = ibObject.second.get<std::string>("geometry.type");
        const Point2D center = Point2D(ibObject.second.get<std::string>("geometry.center"));

        if(type == "circle")
        {
            ibObjs_.push_back(
                        ImmersedBoundaryObject(ibObject.first,
                                               solver_.grid(),
                                               center,
                                               ibObject.second.get<Scalar>("geometry.radius"))
                        );
        }
        else if(type == "box")
        {
            Scalar halfWidth = ibObject.second.get<Scalar>("geometry.width")/2.,
                    halfHeight = ibObject.second.get<Scalar>("geometry.height")/2.;

            std::vector<Point2D> box = {
                center + Vector2D(-halfWidth, -halfHeight),
                center + Vector2D(halfWidth, -halfHeight),
                center + Vector2D(halfWidth, halfHeight),
                center + Vector2D(-halfWidth, halfHeight)
            };

            ibObjs_.push_back(
                        ImmersedBoundaryObject(ibObject.first,
                                               solver_.grid(),
                                               box)
                        );
        }
        else if(type == "polygon")
        {
            std::ifstream fin;
            std::vector<Point2D> verts;
            std::string filename = "case/";
            filename += ibObject.second.get<std::string>("geometry.file").c_str();

            fin.open(filename.c_str());

            if(!fin.is_open())
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "failed to open file \"" + filename + "\".");

            printf("Reading data for \"%s\" from file \"%s\".\n", ibObject.first.c_str(), filename.c_str());

            while(!fin.eof())
            {
                Scalar x, y;
                fin >> x;
                fin >> y;

                verts.push_back(Point2D(x, y));
            }

            fin.close();

            Vector2D translation = center - Polygon(verts).centroid();

            for(Point2D &vert: verts)
                vert += translation;

            ibObjs_.push_back(
                        ImmersedBoundaryObject(ibObject.first,
                                               solver_.grid(),
                                               verts)
                        );
        }
        else
            throw Exception("ImmersedBoundaryObject", "ImmersedBoundaryObject", "invalid geometry type \"" + type + "\".");

        //- Optional geometry parameters
        boost::optional<Scalar> scaleFactor = ibObject.second.get_optional<Scalar>("geometry.scale");
        boost::optional<Scalar> rotationAngle = ibObject.second.get_optional<Scalar>("geometry.rotate");

        if(scaleFactor)
        {
            printf("Scaling \"%s\" by a factor of %lf.\n", ibObject.first.c_str(), scaleFactor.get());
            ibObjs_.back().shape().scale(scaleFactor.get());
        }

        if(rotationAngle)
        {
            printf("Rotating \"%s\" by an angle of %lf degrees.\n", ibObject.first.c_str(), rotationAngle.get());
            ibObjs_.back().shape().rotate(rotationAngle.get()*M_PI/180.);
        }

        //- Boundary information
        for(const auto& child: ibObject.second)
        {
            if(child.first == "geometry")
                continue;

            std::string type = child.second.get<std::string>("type");
            ImmersedBoundaryObject::BoundaryType boundaryType;
            Scalar boundaryRefValue = 0.;

            if(type == "fixed")
                boundaryType = ImmersedBoundaryObject::FIXED;
            else if(type == "normal_gradient")
                boundaryType = ImmersedBoundaryObject::NORMAL_GRADIENT;
            else if(type == "contact_angle")
                boundaryType = ImmersedBoundaryObject::CONTACT_ANGLE;
            else if(type == "partial_slip")
            {
                boundaryType = ImmersedBoundaryObject::PARTIAL_SLIP;
                boundaryRefValue = child.second.get<Scalar>("lambda");
            }
            else
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "unrecognized boundary type \"" + type + "\".");

            printf("Setting boundary type \"%s\" for field \"%s\".\n", type.c_str(), child.first.c_str());

            ibObjs_.back().addBoundaryType(child.first, boundaryType);
            ibObjs_.back().addBoundaryRefValue(child.first, boundaryRefValue);

            if(child.first == "p")
                ibObjs_.back().addBoundaryType("pCorr", boundaryType);
        }

        ibObjs_.back().setInternalCells();
    }

    setCellStatus();
}

Equation<ScalarFiniteVolumeField> ImmersedBoundary::eqns(ScalarFiniteVolumeField &field)
{
    return gc::ib(ibObjs_, field);
}

Equation<VectorFiniteVolumeField> ImmersedBoundary::eqns(VectorFiniteVolumeField &field)
{
    return gc::ib(ibObjs_, field);
}

bool ImmersedBoundary::isIbCell(const Cell &cell) const
{
    for(const ImmersedBoundaryObject& ibObj: ibObjs_)
    {
        if(ibObj.cells().isInGroup(cell))
            return true;
    }

    return false;
}

//- Protected

void ImmersedBoundary::setCellStatus()
{
    for(const Cell &cell: solver_.grid().inactiveCells())
        cellStatus_[cell.id()] = SOLID;

    for(const Cell &cell: solver_.grid().fluidCells())
        cellStatus_[cell.id()] = FLUID;

    for(const ImmersedBoundaryObject& ibObj: ibObjs_)
        for(const Cell &cell: ibObj.cells())
            cellStatus_[cell.id()] = IB;
}
